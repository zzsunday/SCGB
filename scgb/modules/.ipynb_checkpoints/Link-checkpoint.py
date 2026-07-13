#!/usr/bin/env python3
############################################################################
# -*- coding:utf8 -*-
# script: Link.py
# Link plasmid contigs to dereplicated genomes with single-cell reads.
############################################################################

from __future__ import annotations

import csv
import json
import shlex
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Union

import click
import pandas as pd
from Bio import SeqIO
from loguru import logger

from modules.SummaryBc2Cbin import Bc2Ref, GenerateBc2refFile, Summary


def q(value: Union[str, Path]) -> str:
    """Shell-quote a path or scalar."""
    return shlex.quote(str(value))


def run(cmd: str, dry_run: bool = False) -> int:
    """Run a shell command with logging."""
    logger.info("CMD: {}", cmd)
    if dry_run:
        logger.info("Dry run enabled; command was not executed.")
        return 0
    try:
        subprocess.run(cmd, shell=True, check=True)
        return 0
    except subprocess.CalledProcessError as exc:
        logger.error("Command failed with exit code {}: {}", exc.returncode, cmd)
        sys.exit(exc.returncode)


def ensure_file(path: Union[str, Path], desc: str = "file") -> Path:
    p = Path(path).expanduser()
    if not p.is_file():
        raise FileNotFoundError(f"{desc} not found: {p}")
    return p


def ensure_dir(path: Union[str, Path], desc: str = "directory") -> Path:
    p = Path(path).expanduser()
    if not p.is_dir():
        raise NotADirectoryError(f"{desc} not found: {p}")
    return p


def strip_fastq_suffix(path: Path) -> str:
    name = path.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def discover_fastas(fasta_dir: Path, suffixes: tuple[str, ...]) -> list[Path]:
    files: list[Path] = []
    for suffix in suffixes:
        files.extend(sorted(fasta_dir.glob(f"*{suffix}")))
    unique = sorted(dict.fromkeys(files))
    if not unique:
        raise FileNotFoundError(f"No FASTA files found in {fasta_dir} with suffixes {suffixes}")
    return unique


def write_path_list(paths: Iterable[Path], out_path: Path) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fh:
        for path in paths:
            fh.write(f"{path}\n")
    logger.info("Wrote FASTA path list: {}", out_path)
    return out_path


def concatenate_fastas(fasta_paths: Iterable[Path], out_fasta: Path, dry_run: bool = False) -> Path:
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    if dry_run:
        logger.info("Dry run enabled; would concatenate FASTA files into {}", out_fasta)
        return out_fasta
    with out_fasta.open("wb") as fout:
        for fasta in fasta_paths:
            with fasta.open("rb") as fin:
                fout.write(fin.read())
            fout.write(b"\n")
    logger.info("Concatenated dereplicated genomes: {}", out_fasta)
    return out_fasta


def bwa_index(reference: Path, bwa: str, dry_run: bool = False) -> None:
    expected = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    if all((reference.with_suffix(reference.suffix + suffix)).exists() for suffix in expected):
        logger.info("BWA index already exists for {}", reference)
        return
    run(f"{q(bwa)} index {q(reference)}", dry_run=dry_run)


def map_fastq_to_reference(
    reference: Path,
    fastq: Path,
    out_bam: Path,
    threads: int,
    bwa: str,
    samtools: str,
    dry_run: bool = False,
) -> Path:
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    if out_bam.is_file():
        logger.info("Skip genome mapping; BAM already exists: {}", out_bam)
        return out_bam
    bwa_index(reference, bwa=bwa, dry_run=dry_run)
    cmd = (
        f"{q(bwa)} mem -t {threads} {q(reference)} {q(fastq)} "
        f"| {q(samtools)} view -b -F 2308 -q 20 "
        f"| {q(samtools)} sort -@ {threads} -o {q(out_bam)}"
    )
    run(cmd, dry_run=dry_run)
    return out_bam


def normalize_counts_table(counts_file: Path, outdir: Path) -> Path:
    df = pd.read_csv(counts_file, sep=None, engine="python")
    rename = {
        "cell": "CellID",
        "barcode": "CellID",
        "Barcode": "CellID",
        "umi": "Uminum",
        "UMI": "Uminum",
        "UmiNum": "Uminum",
        "reads": "ReadNum",
        "Reads": "ReadNum",
        "readnum": "ReadNum",
    }
    df = df.rename(columns={c: rename.get(c, c) for c in df.columns})
    required = {"CellID", "Uminum", "ReadNum"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"counts file missing columns {missing}; found {list(df.columns)}")
    out_path = outdir / "counts.normalized.tsv"
    df[list(required)].to_csv(out_path, sep="\t", index=False)
    return out_path


def normalize_barcode_list(barcodes_file: Path, outdir: Path) -> tuple[Path, dict[str, str]]:
    """
    Write two-column read-prefix/barcode TSV for Summary and return Barcode -> read prefix.

    Input is expected as two tab-separated columns:
      read_prefix<TAB>barcode
    """
    raw = pd.read_csv(barcodes_file, sep="\t", header=None, dtype=str)
    raw = raw.dropna(how="all")
    out_path = outdir / "barcodes.normalized.tsv"

    if raw.shape[1] >= 2:
        raw = raw.iloc[:, :2].copy()
        raw.columns = ["ReadPrefix", "Barcode"]
        header_mask = (
            raw["ReadPrefix"].str.lower().isin({"index", "idx", "read", "reads", "readprefix", "file"}) |
            raw["Barcode"].str.lower().isin({"barcode", "barcodes", "cellid", "cell"})
        )
        raw = raw.loc[~header_mask].copy()
        raw[["ReadPrefix", "Barcode"]].to_csv(out_path, sep="\t", index=False, header=False)
        prefix_map = dict(zip(raw["Barcode"].astype(str), raw["ReadPrefix"].astype(str)))
    else:
        raise ValueError("barcodes file must be two tab-separated columns: read_prefix<TAB>barcode")

    logger.info("Normalized barcode list: {}", out_path)
    return out_path, prefix_map


def load_meta_fastq_map(meta_file: Path | None) -> dict[str, str]:
    """Load read names from the first column of a tab-separated meta file."""
    if meta_file is None:
        return {}
    meta_file = ensure_file(meta_file, "meta file")
    mapping: dict[str, str] = {}
    with meta_file.open("r", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue
            read_name = row[0].strip()
            if not read_name or read_name.lower() in {"barcode", "prefix", "cell", "read", "reads", "readprefix", "file"}:
                continue
            mapping[read_name] = read_name
    logger.info("Loaded {} read names from meta file", len(mapping))
    return mapping


def candidate_fastqs(single_reads_dir: Path, read_name: str, suffix: str) -> list[Path]:
    """Return candidate FASTQ paths for a read name."""
    candidates = []
    base = Path(read_name).expanduser()
    if base.suffix in {".fq", ".fastq"} or read_name.endswith((".fq.gz", ".fastq.gz")):
        candidates.append(base if base.is_absolute() else single_reads_dir / base)
    else:
        for ext in (suffix, ".fq", ".fastq", ".fq.gz", ".fastq.gz"):
            candidates.append(single_reads_dir / f"{read_name}{ext}")
    return candidates


def barcode_to_fastq(
    barcode: str,
    read_prefix: str,
    single_reads_dir: Path,
    meta_map: dict[str, str],
    suffix: str,
) -> Path:
    for key in (barcode, read_prefix):
        if key in meta_map:
            for fq in candidate_fastqs(single_reads_dir, meta_map[key], suffix):
                if fq.is_file():
                    return fq
            return candidate_fastqs(single_reads_dir, meta_map[key], suffix)[0]
    for fq in candidate_fastqs(single_reads_dir, read_prefix, suffix):
        if fq.is_file():
            return fq
    return candidate_fastqs(single_reads_dir, read_prefix, suffix)[0]


def make_bin_barcode_files(
    summary_bin_dir: Path,
    out_bin_dir: Path,
    barcode_prefix_map: dict[str, str],
) -> list[Path]:
    out_bin_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []

    for source in sorted(summary_bin_dir.glob("filtered_barcodes_*.txt")):
        ref = source.name.removeprefix("filtered_barcodes_").removesuffix(".txt")
        out_file = out_bin_dir / f"{ref}_bc.txt"
        with source.open("r", encoding="utf-8") as fin, out_file.open("w", encoding="utf-8") as fout:
            for line in fin:
                bc = line.strip()
                if not bc:
                    continue
                fout.write(f"{barcode_prefix_map.get(bc, bc)}\n")
        outputs.append(out_file)
        logger.info("Wrote barcode-read-prefix list: {}", out_file)

    if not outputs:
        logger.warning("No filtered_barcodes_*.txt files found in {}", summary_bin_dir)
    return outputs


def concatenate_genome_reads(
    bc2ref_csv: Path,
    outdir: Path,
    single_reads_dir: Path,
    meta_map: dict[str, str],
    barcode_prefix_map: dict[str, str],
    single_read_suffix: str,
    dry_run: bool = False,
) -> Path:
    bc2ref = pd.read_csv(bc2ref_csv, sep="\t")
    required = {"Barcode", "Reference"}
    missing = required - set(bc2ref.columns)
    if missing:
        raise ValueError(f"bc2ref missing columns {missing}; found {list(bc2ref.columns)}")

    reads_dir = outdir / "Reads"
    bin_bc_dir = outdir / "BinBcTxt"
    reads_dir.mkdir(parents=True, exist_ok=True)
    bin_bc_dir.mkdir(parents=True, exist_ok=True)

    for ref, group in bc2ref.groupby("Reference"):
        ref_name = Path(str(ref)).stem
        out_fastq = reads_dir / f"{ref_name}.fq"
        bc_file = bin_bc_dir / f"{ref_name}_bc.txt"
        with bc_file.open("w", encoding="utf-8") as bcf:
            for bc in group["Barcode"].astype(str):
                bcf.write(f"{barcode_prefix_map.get(bc, bc)}\n")

        logger.info("Concatenate reads for {} -> {}", ref, out_fastq)
        if dry_run:
            logger.info("Dry run enabled; would write {}", out_fastq)
            continue

        with out_fastq.open("wb") as fout:
            for bc in group["Barcode"].astype(str):
                read_prefix = barcode_prefix_map.get(bc, bc)
                fq = barcode_to_fastq(
                    barcode=bc,
                    read_prefix=read_prefix,
                    single_reads_dir=single_reads_dir,
                    meta_map=meta_map,
                    suffix=single_read_suffix,
                )
                if not fq.is_file():
                    logger.warning("Missing single-cell FASTQ for barcode {}: {}", bc, fq)
                    continue
                with fq.open("rb") as fin:
                    fout.write(fin.read())
                fout.write(b"\n")

    return reads_dir


def map_genome_reads_to_plasmids(
    reads_dir: Path,
    plasmid_fasta: Path,
    outdir: Path,
    sample: str,
    threads: int,
    bwa: str,
    samtools: str,
    dry_run: bool = False,
) -> Path:
    bam_dir = outdir / "all_plasmid" / "BinReads2Plasimid"
    bam_dir.mkdir(parents=True, exist_ok=True)
    bwa_index(plasmid_fasta, bwa=bwa, dry_run=dry_run)

    for fastq in sorted(reads_dir.glob("*.fq")):
        ref_name = strip_fastq_suffix(fastq)
        bam_out = bam_dir / f"{sample}_p_{ref_name}.bam"
        cmd = (
            f"{q(bwa)} mem -t {threads} -M {q(plasmid_fasta)} {q(fastq)} "
            f"| {q(samtools)} view -b -F 2308 -q 20 "
            f"| {q(samtools)} sort -@ {threads} -o {q(bam_out)}"
        )
        run(cmd, dry_run=dry_run)
    return bam_dir


def bam_to_contact_json(bam_dir: Path, outdir: Path, sample: str, dry_run: bool = False) -> Path:
    try:
        import pysam  # type: ignore
    except Exception:
        logger.error("pysam is required to parse BAM. Install via: pip install pysam")
        sys.exit(1)

    json_dir = outdir / "all_plasmid" / "contact_json"
    json_dir.mkdir(parents=True, exist_ok=True)

    for bam_file in sorted(bam_dir.glob("*.bam")):
        stem = bam_file.stem
        prefix = stem
        if prefix.startswith(f"{sample}_p_"):
            prefix = prefix[len(f"{sample}_p_") :]
        out_json = json_dir / f"{prefix}_p.json"

        if dry_run:
            logger.info("Dry run enabled; would parse {} -> {}", bam_file, out_json)
            continue

        if bam_file.stat().st_size <= 600:
            logger.info("{} exists but is too small (<600B)", bam_file)
            continue

        bin2bc: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for align in bam:
                if align.is_unmapped:
                    continue
                qname = align.query_name or ""
                bc = qname.split("_", 1)[0] if "_" in qname else qname
                contig = align.reference_name
                if contig:
                    bin2bc[contig][bc] += 1

        with out_json.open("w", encoding="utf-8") as fh:
            json.dump({k: dict(v) for k, v in bin2bc.items()}, fh, indent=4)
        logger.info("Wrote contact JSON: {}", out_json)

    return json_dir


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def filter_contact(
    data: dict[str, dict[str, int]],
    purity_df: pd.DataFrame,
    barcode_alias_map: dict[str, str] | None = None,
) -> tuple[dict[str, dict[str, int]], dict[str, dict[str, int]], dict[str, dict[str, list[float]]]]:
    drop_bc: dict[str, dict[str, list[float]]] = defaultdict(dict)
    filter_bc: dict[str, dict[str, int]] = defaultdict(dict)
    purity_dict = purity_df.set_index("Barcode")[["Purity", "ReadNum"]].to_dict(orient="index")
    barcode_alias_map = barcode_alias_map or {}

    for contig, barcodes_dict in data.items():
        for bc, count in barcodes_dict.items():
            info = purity_dict.get(bc)
            if not info and bc in barcode_alias_map:
                info = purity_dict.get(barcode_alias_map[bc])
            if not info:
                continue
            purity_value = float(info["Purity"])
            read_num = float(info["ReadNum"])
            exp_reads = max(5.0, read_num * purity_value * 0.005 / 100.0)

            if 80 <= purity_value <= 99:
                if count < exp_reads:
                    drop_bc[contig][bc] = [count, exp_reads, purity_value, read_num]
                else:
                    filter_bc[contig][bc] = count
            elif purity_value >= 100:
                filter_bc[contig][bc] = count
            else:
                drop_bc[contig][bc] = [count, exp_reads, purity_value, read_num]

    return data, filter_bc, drop_bc


def summarize_contacts(
    json_dir: Path,
    purity_path: Path,
    outdir: Path,
    barcode_alias_map: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, Path]:
    purity = pd.read_csv(purity_path, sep="\t")
    required = {"Barcode", "Purity", "ReadNum"}
    missing = required - set(purity.columns)
    if missing:
        raise ValueError(f"purity table missing columns {missing}; found {list(purity.columns)}")

    thresholded: dict[str, dict[str, int]] = defaultdict(dict)

    for json_file in sorted(json_dir.glob("*_p.json")):
        bin_name = json_file.name.removesuffix("_p.json")
        data, filter_bc, _drop_bc = filter_contact(load_json(json_file), purity, barcode_alias_map)
        filtered_total = sum(sum(barcodes.values()) for barcodes in filter_bc.values())
        filtered_by_contig = {contig: sum(barcodes.values()) for contig, barcodes in filter_bc.items()}
        threshold = filtered_total * 0.01

        for contig, value in filtered_by_contig.items():
            if value > threshold and value > 40:
                thresholded[bin_name][contig] = value

    rows = []
    for bin_name, plasmids in thresholded.items():
        for plasmid_id, readnum in plasmids.items():
            rows.append((bin_name, plasmid_id, readnum))

    df = pd.DataFrame(rows, columns=["Bin", "Plasmid", "ReadNum"])
    out_tsv = outdir / "plasmid_to_genome.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    logger.success("Wrote plasmid-to-genome links: {}", out_tsv)

    out_json = outdir / "plasmid_to_genome.json"
    with out_json.open("w", encoding="utf-8") as fh:
        json.dump({k: dict(v) for k, v in thresholded.items()}, fh, indent=4)
    return df, out_tsv


def build_bin_to_plasmids(ref_df: pd.DataFrame) -> dict[str, list[str]]:
    for col in ("Plasmid", "Bin"):
        if col not in ref_df.columns:
            raise ValueError(f"Missing column: {col}")
    ref_df = ref_df.copy()
    ref_df["Bin_norm"] = ref_df["Bin"].astype(str).str.replace(r"\s+", "", regex=True)
    return (
        ref_df.groupby("Bin_norm")["Plasmid"]
        .apply(lambda x: [str(i) for i in x.dropna().tolist()])
        .to_dict()
    )


def make_plasmid_index(fasta_path: Path):
    return SeqIO.index(str(fasta_path), "fasta")


def find_plasmid_record(idx, key: str):
    if key in idx:
        return idx[key]
    k_ws = key.split()[0]
    if k_ws in idx:
        return idx[k_ws]
    k_pipe = key.split("|")[0]
    if k_pipe in idx:
        return idx[k_pipe]
    parts = key.split("|")
    if parts:
        candidate = parts[-1].split()[0]
        if candidate in idx:
            return idx[candidate]
    return None


def write_genome_with_plasmids(genome_path: Path, records_to_append: list[Any], out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_fp = out_dir / genome_path.name
    with out_fp.open("wb") as fout, genome_path.open("rb") as fin:
        fout.write(fin.read())
        fout.write(b"\n")
    if records_to_append:
        with out_fp.open("a", encoding="utf-8", newline="\n") as fout:
            SeqIO.write(records_to_append, fout, "fasta")
    return out_fp


def append_plasmids_to_genomes(
    link_df: pd.DataFrame,
    genome_dir: Path,
    plasmid_fasta: Path,
    outdir: Path,
    genome_suffix: str,
) -> None:
    bin2plas = build_bin_to_plasmids(link_df)
    plasmid_idx = make_plasmid_index(plasmid_fasta)
    genomes = sorted(genome_dir.glob(f"*{genome_suffix}"))
    if not genomes:
        logger.warning("No genomes found under {} with suffix {}", genome_dir, genome_suffix)

    try:
        for genome in genomes:
            keys_to_try = {genome.stem, genome.name}
            matched_key = None
            for key in keys_to_try:
                normalized = key.replace(" ", "")
                if normalized in bin2plas:
                    matched_key = normalized
                    break

            if matched_key is None:
                out_fp = write_genome_with_plasmids(genome, [], outdir)
                logger.info("[NO PLASMID] Wrote genome only -> {}", out_fp.name)
                continue

            records = []
            missing = []
            for plasmid_header in bin2plas[matched_key]:
                record = find_plasmid_record(plasmid_idx, plasmid_header)
                if record is None:
                    missing.append(plasmid_header)
                else:
                    records.append(record)

            if missing:
                logger.warning(
                    "{}: {} plasmid header(s) not found in DB: {}",
                    genome.name, len(missing), missing[:5]
                )
            out_fp = write_genome_with_plasmids(genome, records, outdir)
            logger.info("[WITH PLASMID] Wrote genome + {} plasmid(s) -> {}", len(records), out_fp.name)
    finally:
        plasmid_idx.close()


@click.command(
    help="Use single-cell reads to link plasmid contigs to dereplicated genomes.",
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option("-fq", "--fastq", required=True,
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              help="All concatenated single-cell sequencing reads FASTQ.")
@click.option("--single-reads-dir", required=True,
              type=click.Path(exists=True, file_okay=False, path_type=Path),
              help="Directory containing per-barcode merged FASTQ files.")
@click.option("-c", "--counts", required=True,
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              help="counts.csv/tsv with CellID, Uminum and ReadNum columns.")
@click.option("-fa_dir", "--fa-dir", required=True,
              type=click.Path(exists=True, file_okay=False, path_type=Path),
              help="Dereplicated genome FASTA directory.")
@click.option("-p", "--plasmid", required=True,
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              help="Plasmid FASTA file.")
@click.option("--meta", required=True,
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              help="Two-column TSV: read_prefix<TAB>barcode.")
@click.option("-o", "--outdir", required=True,
              type=click.Path(file_okay=False, path_type=Path),
              help="Output directory.")
@click.option("--prefix", default=None,
              help="Output prefix. Default: FASTQ basename.")
@click.option("-t", "--threads", default=16, show_default=True,
              type=click.IntRange(min=1), help="Threads for BWA/samtools.")
@click.option("--genome-suffix", default=".fa", show_default=True,
              help="Genome FASTA suffix under --fa-dir.")
@click.option("--single-read-suffix", default=".fq", show_default=True,
              help="First FASTQ suffix tried with read names; .fq/.fastq/.fq.gz/.fastq.gz are also tried.")
@click.option("--dry-run", is_flag=True, help="Print external commands without executing them.")
def main(
    fastq: Path,
    single_reads_dir: Path,
    counts: Path,
    fa_dir: Path,
    plasmid: Path,
    meta: Path,
    outdir: Path,
    prefix: str | None,
    threads: int,
    genome_suffix: str,
    single_read_suffix: str,
    dry_run: bool,
) -> None:
    outdir = Path(outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    logger.add(outdir / "link.log",
               level="INFO",
               format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}",
               mode="a",
               enqueue=True)

    fastq = ensure_file(fastq, "single-cell FASTQ")
    single_reads_dir = ensure_dir(single_reads_dir, "single-cell reads directory")
    fa_dir = ensure_dir(fa_dir, "dereplicated FASTA directory")
    plasmid = ensure_file(plasmid, "plasmid FASTA")
    sample_name = prefix or strip_fastq_suffix(fastq)
    sample_dir = outdir / sample_name
    sample_dir.mkdir(parents=True, exist_ok=True)
    bwa = "bwa"
    samtools = "samtools"

    fasta_paths = discover_fastas(fa_dir, (genome_suffix,))
    drep_ref_fasta = outdir / "all_c.fa"
    concatenate_fastas(fasta_paths, drep_ref_fasta, dry_run=dry_run)

    path_list = write_path_list(fasta_paths, outdir / "all_binette_bin_c_path.txt")
    genome2id, fasta_n, contig_n = GenerateBc2refFile(path_list, outdir)
    logger.info("Genome2id: {} | genomes={} contigs={}", genome2id, fasta_n, contig_n)

    genome_bam = sample_dir / f"{sample_name}_c_mapped_sorted.bam"
    map_fastq_to_reference(
        reference=drep_ref_fasta,
        fastq=fastq,
        out_bam=genome_bam,
        threads=threads,
        bwa=bwa,
        samtools=samtools,
        dry_run=dry_run,
    )

    if dry_run:
        logger.info("Dry run stops before Python parsing steps that require BAM outputs.")
        return

    bc2ref_csv, detail_json = Bc2Ref(bam_path=genome_bam, ref_map_path=genome2id, outdir=sample_dir)
    counts_tsv = normalize_counts_table(counts, sample_dir)
    barcodes_txt, barcode_prefix_map = normalize_barcode_list(meta, sample_dir)
    Summary(
        counts_file=counts_tsv,
        filter_bc_file=barcodes_txt,
        outdir=sample_dir,
        bc2ref_csv=bc2ref_csv,
        detail_json=detail_json,
    )

    make_bin_barcode_files(sample_dir / "BinBcTxt", sample_dir / "BinBcTxt", barcode_prefix_map)
    meta_map = load_meta_fastq_map(meta)
    reads_dir = concatenate_genome_reads(
        bc2ref_csv=Path(bc2ref_csv),
        outdir=sample_dir,
        single_reads_dir=single_reads_dir,
        meta_map=meta_map,
        barcode_prefix_map=barcode_prefix_map,
        single_read_suffix=single_read_suffix,
        dry_run=dry_run,
    )

    plasmid_bam_dir = map_genome_reads_to_plasmids(
        reads_dir=reads_dir,
        plasmid_fasta=plasmid,
        outdir=sample_dir,
        sample=sample_name,
        threads=threads,
        bwa=bwa,
        samtools=samtools,
        dry_run=dry_run,
    )
    contact_json_dir = bam_to_contact_json(plasmid_bam_dir, sample_dir, sample_name, dry_run=dry_run)
    prefix_to_barcode = {prefix: barcode for barcode, prefix in barcode_prefix_map.items()}
    link_df, link_tsv = summarize_contacts(
        contact_json_dir,
        sample_dir / "purity.tsv",
        sample_dir,
        barcode_alias_map=prefix_to_barcode,
    )

    append_plasmids_to_genomes(
        link_df=link_df,
        genome_dir=fa_dir,
        plasmid_fasta=plasmid,
        outdir=outdir / "genomes",
        genome_suffix=genome_suffix,
    )

    logger.success("Link workflow finished. Plasmid links: {}", link_tsv)


if __name__ == "__main__":
    main()
