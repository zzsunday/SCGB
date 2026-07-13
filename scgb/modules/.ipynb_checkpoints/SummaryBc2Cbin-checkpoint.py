#!/usr/bin/env python
############################################################################
# -*- coding:utf-8 -*-
# version: 1.0
# author: sundayzz
# mail: sundayzz@sina.cn
# script:SummaryBc2Cbin.py
# Script: Summary barcode infomation, purity, mean reads...
############################################################################

"""
Script: Summarize barcode information (optionally reference-aware purity).
Pipeline:
  1) Scan FASTA bins -> write <outdir>/bin_paths.txt
  2) Build <outdir>/Genome2id.txt (contig -> fasta basename)
  3) From BAM + Genome2id, compute:
       - <outdir>/tmp.csv        (Barcode, Reference, Count)
       - <outdir>/bc2ref.csv     (best Reference per Barcode)
       - <outdir>/tmp_detail.json
  4) Summarize counts + (optional) purity -> <outdir>/summary.tsv/json
Notes:
  - Code/comments in English per your convention.
"""

from __future__ import annotations
import os
import sys
import json
from pathlib import Path
from typing import Union, Iterable
from collections import defaultdict, Counter
import click
from loguru import logger

# ---------- Optional dependencies ----------
try:
    from Bio import SeqIO  # type: ignore
except Exception:
    SeqIO = None  # validated at runtime

try:
    import pandas as pd  # type: ignore
except Exception:
    pd = None  # validated at runtime

try:
    import pysam  # type: ignore
except Exception:
    pysam = None  # validated at runtime

# ---------- Constants ----------
EXTS = {".fa", ".fna", ".fasta"}
FIXED_OUTNAME = "bin_paths.txt"
GENOME2ID_NAME = "Genome2id.txt"


# ---------- Guards / validators ----------
def ensure_dependency(mod, name: str, install_hint: str) -> None:
    """Ensure a Python module is importable; exit with a helpful hint if missing."""
    if mod is None:
        logger.error("{} is required but not installed. {}", name, install_hint)
        sys.exit(1)


def ensure_dir(path: Union[str, Path]) -> Path:
    """Ensure a directory exists and is writable."""
    p = Path(path).expanduser().resolve()
    p.mkdir(parents=True, exist_ok=True)
    try:
        test_file = p / ".write_test"
        with test_file.open("w") as fh:
            fh.write("ok")
        test_file.unlink(missing_ok=True)
    except Exception as e:
        logger.error("Directory is not writable: {} ({})", p, e)
        sys.exit(1)
    return p


def ensure_file(path: Union[str, Path], desc: str = "file") -> Path:
    """Ensure a readable file exists."""
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        logger.error("{} not found: {}", desc, p)
        sys.exit(1)
    try:
        with p.open("rb"):
            pass
    except Exception as e:
        logger.error("Cannot read {}: {} ({})", desc, p, e)
        sys.exit(1)
    return p


# ---------- FASTA scanning & mapping ----------
def get_bin_info(indir: Union[str, Path], outdir: Union[str, Path], recursive: bool = True) -> tuple[str, int]:
    """Scan a folder for FASTA-like files and write absolute paths to <outdir>/bin_paths.txt."""
    root = Path(indir).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"Input directory does not exist or is not a directory: {root}")

    paths: list[str] = []
    if recursive:
        for dirpath, _, filenames in os.walk(root):
            for fn in filenames:
                if Path(fn).suffix.lower() in EXTS:
                    p = Path(dirpath) / fn
                    paths.append(str(p.resolve()))
    else:
        for p in root.iterdir():
            if p.is_file() and p.suffix.lower() in EXTS:
                paths.append(str(p.resolve()))

    unique_sorted = sorted(set(paths))

    out_dir = ensure_dir(outdir)
    out_file = out_dir / FIXED_OUTNAME

    with out_file.open("w", encoding="utf-8") as fh:
        fh.write("\n".join(unique_sorted))

    return str(out_file), len(unique_sorted)


def GenerateBc2refFile(out_path: Union[str, Path], outdir: Union[str, Path]) -> tuple[str, int, int]:
    """Build '<outdir>/Genome2id.txt' mapping contig -> fasta basename from a list of fasta paths."""
    ensure_dependency(SeqIO, "Biopython", "Install via: pip install biopython")

    list_file = ensure_file(out_path, desc="path list file")
    out_dir = ensure_dir(outdir)
    output = out_dir / GENOME2ID_NAME

    fasta_count = 0
    contig_count = 0

    with list_file.open("r", encoding="utf-8") as fh_in, output.open("w", encoding="utf-8") as fh_out:
        for line in fh_in:
            fpath_str = line.strip()
            if not fpath_str:
                continue
            fpath = Path(fpath_str)
            if not fpath.is_file():
                logger.warning("Listed FASTA does not exist, skip: {}", fpath_str)
                continue

            fasta_count += 1
            prefix = os.path.basename(fpath_str)  # keep extension to avoid accidental collisions

            try:
                with fpath.open("r", encoding="utf-8") as fh_fa:
                    has_any = False
                    for record in SeqIO.parse(fh_fa, "fasta"):
                        has_any = True
                        fh_out.write(f"{record.id}\t{prefix}\n")
                        contig_count += 1
                    if not has_any:
                        logger.warning("FASTA contains no records (empty?): {}", fpath)
            except Exception as e:
                logger.exception("Failed parsing FASTA: {} ({})", fpath_str, e)

    if fasta_count == 0:
        logger.error("No readable FASTA files were processed; Genome2id.txt would be empty.")
        sys.exit(1)

    if contig_count == 0:
        logger.error("No contigs were written to Genome2id.txt; please check your FASTA inputs.")
        sys.exit(1)

    return str(output), fasta_count, contig_count


# ---------- BAM -> bc2ref ----------
def _validate_bam_vs_map(bam: "pysam.AlignmentFile", ref_map_path: Path) -> None:
    """Sanity-check overlap between BAM reference names and Genome2id contig IDs."""
    ensure_dependency(pd, "pandas", "Install via: pip install pandas")
    genome = pd.read_csv(ref_map_path, sep="\t", header=None, index_col=0)
    contigs_in_map = set(genome.index.astype(str).tolist())

    bam_refs = set(map(str, (bam.header.references or [])))
    if not bam_refs:
        logger.warning("No SQ/contig names found in BAM header.")

    overlap = contigs_in_map & bam_refs
    if not overlap:
        logger.error(
            "No overlap between BAM reference names ({} entries) and Genome2id contigs ({} entries). "
            "Please ensure they are built from the same FASTA set.",
            len(bam_refs), len(contigs_in_map)
        )
        sys.exit(1)
    else:
        logger.info("BAM/Genome2id overlap OK: {} shared contigs", len(overlap))


def Bc2Ref(
    bam_path: Union[str, Path],
    ref_map_path: Union[str, Path],
    outdir: Union[str, Path],
    qname_delim: str = "_",
) -> tuple[str, str]:
    """
    Build Barcode->Reference statistics from an alignment BAM and a contig->reference mapping.

    Returns:
      (bc2ref_csv_path, tmp_detail_json_path)
    """
    ensure_dependency(pd, "pandas", "Install via: pip install pandas")
    ensure_dependency(pysam, "pysam", "Install via: pip install pysam")

    bam_path = ensure_file(bam_path, desc="BAM")
    ref_map_path = ensure_file(ref_map_path, desc="Genome2id.txt")
    outdir = ensure_dir(outdir)

    tmp_file = outdir / "tmp.csv"
    bc2ref_file = outdir / "bc2ref.csv"
    json_file = outdir / "tmp_detail.json"

    logger.info("Write raw table to: {}", tmp_file)
    logger.info("Write bc2ref to   : {}", bc2ref_file)

    # contig_id -> reference_label
    genome = pd.read_csv(ref_map_path, sep="\t", header=None, index_col=0)
    if genome.empty:
        logger.error("Genome2id.txt is empty: {}", ref_map_path)
        sys.exit(1)
    genome2id = genome.to_dict()[1]  # {contig: ref_label}

    d = defaultdict(lambda: defaultdict(int))

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        _validate_bam_vs_map(bam, ref_map_path)

        processed = skipped_unmapped = skipped_no_qname = skipped_no_delim = skipped_ref_not_found = 0

        for align in bam:
            processed += 1
            if align.is_unmapped:
                skipped_unmapped += 1
                continue

            qname = align.query_name
            if not qname:
                skipped_no_qname += 1
                continue

            if qname_delim not in qname:
                skipped_no_delim += 1
                continue

            bc = qname.split(qname_delim, 1)[0]
            ref_id = align.reference_name
            if ref_id is None:
                continue

            ref_label = genome2id.get(ref_id)
            if ref_label is None:
                skipped_ref_not_found += 1
                continue

            d[bc][ref_label] += 1

    logger.info(
        "BAM processed: total={}, unmapped={}, no_qname={}, no_delim='{}'={}, ref_not_found={}",
        processed, skipped_unmapped, skipped_no_qname, qname_delim, skipped_no_delim, skipped_ref_not_found
    )

    # long table
    rows = []
    for bc, ref_counts in d.items():
        for ref_label, count in ref_counts.items():
            rows.append((bc, ref_label, count))

    if not rows:
        logger.error(
            "No usable alignments produced any (Barcode, Reference) counts. "
            "Check BAM qname format, delimiter '{}', and Genome2id mapping.",
            qname_delim
        )
        sys.exit(1)

    pd.DataFrame(rows, columns=["Barcode", "Reference", "Count"]).to_csv(tmp_file, sep="\t", index=False)

    # nested json
    with json_file.open("w", encoding="utf-8") as jf:
        json.dump({bc: dict(ref_counts) for bc, ref_counts in d.items()}, jf, indent=4)

    # best reference per barcode
    bc2ref = {bc: max(ref_counts, key=ref_counts.get) for bc, ref_counts in d.items()}
    pd.DataFrame(list(bc2ref.items()), columns=["Barcode", "Reference"]).to_csv(bc2ref_file, sep="\t", index=False)

    logger.info("bc2ref entries: {}", len(bc2ref))
    logger.info("Done.")
    return str(bc2ref_file), str(json_file)


# ---------- Summary ----------
def Summary(counts_file: Union[str, Path],
            filter_bc_file: Union[str, Path],
            outdir: Union[str, Path],
            bc2ref_csv: Union[str, Path] | None = None,
            detail_json: Union[str, Path] | None = None) -> None:
    """
    Summarize per-cell metrics and (optionally) purity by reference.

    Inputs:
      counts_file  : TSV with columns: CellID, Uminum, ReadNum
      filter_bc_file: TSV with read-prefix and barcode columns
      bc2ref_csv   : optional TSV (Barcode, Reference)
      detail_json  : optional JSON with nested counts per barcode {ref: count}
    """
    outdir = ensure_dir(outdir)
    counts_file = ensure_file(counts_file, "counts file")
    filter_bc_file = ensure_file(filter_bc_file, "filtered barcodes file")

    output_json = Path(outdir) / "summary.json"
    output_tsv = Path(outdir) / "summary.tsv"

    # Load tables
    counts = pd.read_csv(counts_file, sep="\t")
    required_cols = {"CellID", "Uminum", "ReadNum"}
    if not required_cols.issubset(counts.columns):
        logger.error("counts file must contain columns: {}", ", ".join(sorted(required_cols)))
        sys.exit(1)

    filter_bc = pd.read_csv(filter_bc_file, sep="\t", header=None, dtype=str)
    if filter_bc.shape[1] < 2:
        logger.error("filtered barcodes file must have two tab-separated columns: read_prefix and Barcode.")
        sys.exit(1)
    filter_bc = filter_bc.iloc[:, :2].copy()
    filter_bc.columns = ["ReadPrefix", "Barcode"]
    header_mask = (
        filter_bc["ReadPrefix"].str.lower().isin({"index", "idx", "read", "reads", "readprefix", "file"}) |
        filter_bc["Barcode"].str.lower().isin({"barcode", "barcodes", "cellid", "cell"})
    )
    filter_bc = filter_bc.loc[~header_mask].copy()

    # Subset counts to barcodes of interest
    filter_bc_counts = counts[counts["CellID"].isin(filter_bc["Barcode"])].copy()

    # Default: reference-free summary
    filter_cal_ = filter_bc_counts.rename(columns={"CellID": "Barcode"})

    # Optional: reference-aware purity
    if bc2ref_csv and detail_json:
        bc2ref_csv = ensure_file(bc2ref_csv, "bc2ref csv")
        detail_json = ensure_file(detail_json, "detail json")

        bc2ref = pd.read_csv(bc2ref_csv, sep="\t")
        if not {"Barcode", "Reference"}.issubset(bc2ref.columns):
            logger.error("bc2ref csv must contain columns: Barcode, Reference")
            sys.exit(1)

        # purity per barcode from nested JSON
        with open(detail_json, "r") as jf:
            nested = json.load(jf)

        purity_map = {}
        for bc in filter_bc["Barcode"].values:
            counts_dict = nested.get(str(bc)) or nested.get(bc)
            if not counts_dict:
                continue
            total = sum(counts_dict.values())
            if total <= 0:
                continue
            purity = round(max(counts_dict.values()) / total * 100, 0)
            purity_map[bc] = purity

        purity_df = pd.DataFrame(list(purity_map.items()), columns=["Barcode", "Purity"])
        filter_cal_ = purity_df.merge(filter_bc_counts, left_on="Barcode", right_on="CellID", how="inner")
        filter_cal_ = filter_cal_.drop(columns=["CellID"]).merge(bc2ref, on="Barcode", how="inner")


        # write per-reference barcode lists
        BinBcTxt_dir = Path(outdir)/"BinBcTxt"
        BinBcTxt_dir.mkdir(parents=True, exist_ok=True)
        for ref_name in sorted(filter_cal_["Reference"].unique()):
            out_path = BinBcTxt_dir / f"filtered_barcodes_{ref_name}.txt"
            filter_cal_.loc[filter_cal_["Reference"] == ref_name, "Barcode"].to_csv(out_path, index=False, header=False)

        # also write purity table
        (Path(outdir) / "purity.tsv").write_text(
            filter_cal_.to_csv(sep="\t", index=False), encoding="utf-8"
        )

    # Summary statistics (safe guards)
    n_cells = int(filter_cal_.shape[0])
    reads_sum = float(filter_cal_["ReadNum"].sum()) if n_cells else 0.0
    umis_sum = float(filter_cal_["Uminum"].sum()) if n_cells else 0.0
    reads_all = float(counts["ReadNum"].sum()) if counts.shape[0] else 0.0

    summary_dict = {}
    summary_dict["Estimated Barcodes"] = n_cells
    summary_dict["Mean Reads per Cell"] = int(round(reads_sum / n_cells)) if n_cells else 0
    summary_dict["Mean Umis per Cell"] = int(round(umis_sum / n_cells)) if n_cells else 0
    summary_dict["Median Umis per Cell"] = int(round(float(filter_cal_["Uminum"].median()))) if n_cells else 0
    summary_dict["Fraction Reads in Cells (%)"] = int(round((reads_sum / reads_all * 100))) if reads_all > 0 else 0
    summary_dict["Sequencing Saturation (%)"] = int(round((umis_sum / reads_sum * 100))) if reads_sum > 0 else 0

    # Reference-aware extras
    if "Reference" in filter_cal_.columns:
        ref_counts = Counter(filter_cal_["Reference"])
        summary_dict["Reference Counts"] = dict(ref_counts)
        if "Purity" in filter_cal_.columns and n_cells:
            summary_dict["Median purity (%)"] = int(round(float(filter_cal_["Purity"].median())))
            summary_dict["Mean purity (%)"] = int(round(float(filter_cal_["Purity"].mean())))
            summary_dict["Purity over 80%"] = int((filter_cal_["Purity"] >= 80).sum())
            summary_dict["Purity over 90%"] = int((filter_cal_["Purity"] >= 90).sum())
            summary_dict["Purity over 95%"] = int((filter_cal_["Purity"] >= 95).sum())

    # Write outputs
    with open(output_json, "w", encoding="utf-8") as jf:
        json.dump(summary_dict, jf, indent=4)
    pd.DataFrame(list(summary_dict.items()), columns=["Metric", "Value"]).to_csv(output_tsv, sep="\t", index=False)

    logger.info("Summary written: {}", output_tsv)
    logger.info("Summary JSON   : {}", output_json)


# ---------- CLI ----------
@click.command(
    help="Scan FASTA files -> <outdir>/bin_paths.txt, then build <outdir>/Genome2id.txt, "
         "then run BAM->bc2ref, finally compute summary (with or without reference-aware purity).",
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option(
    "-i", "--indir",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
    help="Input folder containing .fa/.fna/.fasta files.",
)
@click.option(
    "-b", "--bam",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, path_type=Path),
    help="Aligned BAM file (query name should encode barcode as 'BARCODE_...').",
)
@click.option(
    "-C", "--counts",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, path_type=Path),
    help="Counts TSV with columns: CellID, Uminum, ReadNum.",
)
@click.option(
    "-B", "--barcode",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, path_type=Path),
    help="Filtered barcodes file (one column, no header).",
)
@click.option(
    "-o", "--outdir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False, dir_okay=True, path_type=Path),
    help="Output directory.",
)
@click.option("--no-recursive", is_flag=True, default=False, help="Only scan the top-level directory (do not recurse).")
@click.option("--qname-delim", default="_", show_default=True,
              help="Delimiter in BAM query name that separates BARCODE from the rest.")
def main(indir: Path, bam: Path, counts: Path, barcode: Path, outdir: Path, no_recursive: bool, qname_delim: str) -> None:
    # logging
    outdir = ensure_dir(outdir)
    logger.remove()
    logger.add(sys.stdout, level="INFO",
               format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level}</level> | {message}")
    logger.add(outdir / "pipeline.log", level="INFO",
               format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}", mode="a", enqueue=True)

    # 1) scan & write bin_paths.txt
    logger.info("Start scanning: indir='{}', recursive={}", indir, not no_recursive)
    out_path, count = get_bin_info(indir=indir, outdir=outdir, recursive=not no_recursive)
    if count == 0:
        logger.error("No .fa/.fna/.fasta files found under: {}", indir)
        sys.exit(1)
    logger.info("Wrote {} path(s) to: {}", count, out_path)

    # 2) build Genome2id.txt
    ensure_dependency(SeqIO, "Biopython", "Install via: pip install biopython")
    gmap_path, fasta_n, contig_n = GenerateBc2refFile(out_path=out_path, outdir=outdir)
    logger.info("Genome2id generated: {}", gmap_path)
    logger.info("Processed FASTA files: {}; Total contigs written: {}", fasta_n, contig_n)

    # 3) bam -> bc2ref
    ensure_dependency(pd, "pandas", "Install via: pip install pandas")
    ensure_dependency(pysam, "pysam", "Install via: pip install pysam")
    bc2ref_csv, detail_json = Bc2Ref(bam_path=bam, ref_map_path=gmap_path, outdir=outdir, qname_delim=qname_delim)

    # 4) summary (reference-aware because we have bc2ref_csv & detail_json)
    Summary(counts_file=counts, filter_bc_file=barcode, outdir=outdir,
            bc2ref_csv=bc2ref_csv, detail_json=detail_json)

    logger.success("All steps finished successfully.")


if __name__ == "__main__":
    main()
