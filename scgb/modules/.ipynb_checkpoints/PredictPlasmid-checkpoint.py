############################################################################
# -*- coding:utf8 -*-
# data: 2025-08-10
# author: sundayzz
# mail: sundayzz@sina.cn
# this is part of SingleCellMeta assembler
# this script is mainly for predicting contig for assembly contigs 
############################################################################

import sys
import shlex
import subprocess
from pathlib import Path
from typing import Union, Tuple
import click
from loguru import logger
from Bio import SeqIO
import pandas as pd

from scgb.modules.FilterFalsePlasmid import split_plasmid_and_chromosome

# =============== check tools ===============
def ensure_file(path_str: Union[str, Path]) -> str:
    """check file exist"""
    p = Path(path_str).expanduser()
    if p.is_file():
        return str(p)
    raise FileNotFoundError(f"File not found: {p}")


def q(value: Union[str, Path]) -> str:
    """Shell-quote a path or scalar."""
    return shlex.quote(str(value))


def tool_cmd(command: Union[str, Path], env_manager: str, env_name: str | None) -> str:
    """
    Build a command that optionally runs inside a conda/micromamba environment.

    If env_name is empty, the command is executed directly from PATH.
    """
    command_str = str(command)
    if not env_name:
        return q(command_str)
    return f"{q(env_manager)} run -n {q(env_name)} {q(command_str)}"


def run(cmd: str, dry_run: bool = False) -> int:
    """run cmd"""
    logger.info("CMD: {}", cmd)
    if dry_run:
        logger.info("Dry run enabled; command was not executed.")
        return 0
    try:
        subprocess.run(cmd, shell=True, check=True)
        return 0
    except subprocess.CalledProcessError as e:
        logger.error("Command failed with exit code {}: {}", e.returncode, cmd)
        sys.exit(e.returncode)


# =============== core function ===============
def filter_contigs_by_len(
    input_file: Union[str, Path],
    k: int,
    write_length_table: bool = True,
    outdir: Union[str, Path] = None,
) -> Path:
    """
    Write contigs to FASTA sequentially using Biopython and keep only records with length >= k.

    - If k == 0, no filtering is performed and the original input file path is returned.
    - The output FASTA is named '<stem>_<k>.fa' and is written to 'outdir' (defaults to the input's directory).
    - If write_length_table is True, also write 'contig_length_filter<k>.txt' with lines 'id<TAB>length'.

    Returns
    -------
    Path
        Path to the output FASTA file.
    """
    in_path = Path(input_file).expanduser()
    if not in_path.is_file():
        raise FileNotFoundError(f"Input FASTA not found: {in_path}")

    if k <= 0:
        logger.info("len==0, skip filtering; use original FASTA: {}", in_path)
        return in_path

    out_base = Path(outdir).expanduser() if outdir else in_path.parent
    out_base.mkdir(parents=True, exist_ok=True)

    out_fa = out_base / f"{in_path.stem}_{k}.fa"
    out_len = out_base / f"contig_length_filter{k}.txt"

    kept = total = 0
    with in_path.open("r") as fin, out_fa.open("w") as fout_fa:
        fout_len = out_len.open("w") if write_length_table else None
        try:
            for rec in SeqIO.parse(fin, "fasta"):
                total += 1
                L = len(rec)
                if L >= k:
                    SeqIO.write(rec, fout_fa, "fasta")  
                    if fout_len:
                        fout_len.write(f"{rec.id}\t{L}\n")
                    kept += 1
        finally:
            if fout_len:
                fout_len.close()

    logger.info("Filtering done: kept {}/{} contigs >= {} bp -> {}", kept, total, k, out_fa)
    return out_fa


def build_genomad_cmd(
    genomad: Union[str, Path],
    input_fa: Union[str, Path],
    db_dir: Union[str, Path],
    threads: int,
    outdir: Union[str, Path],
    cleanup: bool = False,
    extra: Tuple[str, ...] = (),
) -> str:
    """
    construt geNomad cmd：
      genomad end-to-end --threads N  -v <input_fasta> <output_dir> <database_dir> 
    """
    parts = [
        str(genomad),
        "end-to-end",
        "-t", str(int(threads)),
        "-v", q(Path(input_fa).expanduser()),
        q(Path(outdir).expanduser()),
        q(Path(db_dir).expanduser()),
    ]
    if cleanup:
        parts.append("--cleanup")
    if extra:
        parts.extend(extra)
    return " ".join(parts)


def build_plasmer_cmd(
    plasmer: Union[str, Path],
    input_fa: Union[str, Path],
    db_dir: Union[str, Path],
    threads: int,
    outdir: Union[str, Path],
    prefix: str
) -> str:
    """
    p：
      plasmer -v <input_fasta> -o <output_dir> -d <database_dir> -t N -p <prefx>
    """
    parts = [
        str(plasmer),
        "-g", q(Path(input_fa).expanduser()),
        "-t", str(int(threads)),
        "-p", q(prefix),
        "-o", q(Path(outdir).expanduser()),
        "-d", q(Path(db_dir).expanduser()),

    ]
    return " ".join(parts)

def Predict(
    contig_fa: Union[str, Path],
    genomad_path: Union[str, Path],
    plasmer_path: Union[str, Path],
    db1: Union[str, Path],
    db2: Union[str, Path],
    threads: int,
    outdir: Union[str, Path],
    prefix: str,
    env_manager: str,
    genomad_env: str,
    plasmer_env: str,
    dry_run: bool = False,
) -> Tuple[Path, Path]:
    """
    run geNomad and Plasmer，return output path (genomad_out, plasmer_out)。
    """
    outdir = Path(outdir).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)

    genomad_out = outdir / "genomad_out"
    plasmer_out = outdir / "plasmer_out"
    genomad_out.mkdir(parents=True, exist_ok=True)
    plasmer_out.mkdir(parents=True, exist_ok=True)

    genomad_cmd = tool_cmd(genomad_path, env_manager, genomad_env)
    plasmer_cmd = tool_cmd(plasmer_path, env_manager, plasmer_env)

    # geNomad
    cmd1 = build_genomad_cmd(genomad_cmd, contig_fa, db1, threads, genomad_out)
    logger.info(f" genomad cmd: {cmd1}")
    run(cmd1, dry_run=dry_run)

    # Plasmer
    cmd2 = build_plasmer_cmd(plasmer_cmd, contig_fa, db2, threads, plasmer_out, prefix)
    logger.info(f" plasmer cmd: {cmd2}")
    run(cmd2, dry_run=dry_run)

    return genomad_out, plasmer_out


def ExtractPlasmidContigName(outdir, genomad_out, plasmer_out) -> Path:
    
    ## get plasmid contig name from genomad output ##
    matches = sorted(genomad_out.rglob("*plasmid_summary.tsv"))
    if not matches:
        raise FileNotFoundError(f"No '*plasmid_summary.tsv' under: {genomad_out}")
    if len(matches) > 1:
        logger.warning("Multiple matches found: {}", matches)
    genomad_file = matches[0]
    logger.info(f"genomad file: {genomad_file}")
    


    ## get plasmid contig name from plasmer output ##
    target_dir = plasmer_out/ "results" 
    matches = sorted(target_dir.glob("*plasmer.predClass.tsv"))
    if not matches:
        raise FileNotFoundError(f"No '*plasmer.predClass.tsv' under: {target_dir}")
    if len(matches) > 1:
        logger.warning("Multiple matches found: {}", matches)
    plasmer_file = matches[0]
    logger.info(f"plasmer file: {plasmer_file}")
    
   
    dat = pd.read_csv(genomad_file, sep="\t")
    dat1 = pd.read_csv(plasmer_file, sep="\t", header=None)
    plasmid_names = outdir / "plasmid_contig.txt"
    with open(plasmid_names, "w+") as f:
        for i in set(dat["seq_name"]) | set(dat1[dat1[1] == "plasmid"][0]):
            f.write(i+"\n")
    logger.info("Plasmid contig name list: {}", plasmid_names)
    return plasmid_names

# =============== run function ===============
@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option("--infile", "-i",
              required=True,
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              help="Assembly contig FASTA")
@click.option("--min_len", "-l", "min_len",
              default=5000, show_default=True,
              type=click.IntRange(min=0),
              help="Filter contig length (bp). 0 = no filter")
@click.option("--p1", "genomad_path",
              default="genomad", show_default=True,
              help="geNomad executable or command name.")
@click.option("--p2", "plasmer_path",
              default="Plasmer", show_default=True,
              help="Plasmer executable or command name.")
@click.option("--db1", "genomad_db",
              required=True, 
              type=click.Path(exists=True, path_type=Path),
              help="geNomad database directory")
@click.option("--db2", "plasmer_db",
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help="Plasmer database directory")
@click.option("--thread", "-t",
              default=16, show_default=True,
              type=click.IntRange(min=1),
              help="Threads for predictors")
@click.option("--outdir", "-o",
              required=True,
              type=click.Path(file_okay=False, path_type=Path),
              help="Output directory")
@click.option("--prefix", "-p", required=True, help="output file prefix")
@click.option("--seqkit", default="seqkit", show_default=True,
              help="seqkit executable used to split plasmid/chromosome FASTA files.")
@click.option("--env-manager", default="micromamba", show_default=True,
              help="Environment manager used to run geNomad and Plasmer.")
@click.option("--genomad-env", default="scgb_genomad", show_default=True,
              help="Environment containing geNomad. Empty string = use active PATH.")
@click.option("--plasmer-env", default="scgb_plasmer", show_default=True,
              help="Environment containing Plasmer. Empty string = use active PATH.")
@click.option(
    "--extract",
    "-e",
    is_flag=True,
    help="Deprecated compatibility flag; contig splitting is now always performed unless --dry-run.",
)
@click.option("--dry-run", is_flag=True, help="Print geNomad/Plasmer commands without executing them.")

def main(
    infile,
    min_len,
    genomad_path,
    plasmer_path,
    genomad_db,
    plasmer_db,
    thread,
    outdir,
    prefix,
    seqkit,
    env_manager,
    genomad_env,
    plasmer_env,
    extract,
    dry_run,
):
    """
    Predict plasmid-like contigs from a metagenomic assembly using geNomad and Plasmer.

    This is the first simplified SCGB module. It writes:
      - plasmid_contig.txt
      - <prefix>.plasmid_contigs.fa
      - <prefix>.chromosome_contigs.fa

    Workflow:
      1) (Optional) filter user-supplied contigs by length;
      2) Run geNomad & Plasmer;
      3) Split filtered assembly contigs into plasmid and chromosome FASTA files.
    """
    # log file output dir
    outdir.mkdir(parents=True, exist_ok=True)
    logger.remove()
    logger.add(sys.stdout, level="INFO",
               format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level}</level> | {message}")
    logger.add(outdir / "predict.log", level="INFO",
               format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}",
               mode="a", enqueue=True)

    logger.info("Input FASTA: {}", infile)
    logger.info("Min length (bp): {}", min_len)
    logger.info("geNomad: {} | DB: {}", genomad_path, genomad_db)
    logger.info("Plasmer: {} | DB: {}", plasmer_path, plasmer_db)
    logger.info("Outdir: {} | Threads: {}", outdir, thread)
    logger.info("Environment manager: {}", env_manager)
    logger.info("Predictor environments: genomad={} plasmer={}",
                genomad_env or "PATH", plasmer_env or "PATH")

    # 1) filter contig length (optional)
    filtered_fa = infile
    if min_len > 0 :
        filtered_fa = filter_contigs_by_len(infile, min_len, write_length_table=True, outdir=outdir)
        logger.info("FASTA for prediction: {}", filtered_fa)

    # 2) predict
    genomad_out, plasmer_out = Predict(
        contig_fa=filtered_fa,
        genomad_path=genomad_path,
        plasmer_path=plasmer_path,
        db1=genomad_db,
        db2=plasmer_db,
        threads=thread,
        outdir=outdir,
        prefix=prefix,
        env_manager=env_manager,
        genomad_env=genomad_env,
        plasmer_env=plasmer_env,
        dry_run=dry_run

    )

    logger.info("geNomad out: {}", genomad_out)
    logger.info("Plasmer out: {}", plasmer_out)   
    if not dry_run:
        logger.info("===== Extract plasmid contigs and split FASTA outputs =====")
        plasmid_names = ExtractPlasmidContigName(outdir, genomad_out, plasmer_out)
        plasmid_fa = outdir / f"{prefix}.plasmid_contigs.fa"
        chrom_fa = outdir / f"{prefix}.chromosome_contigs.fa"
        split_plasmid_and_chromosome(
            names_file=plasmid_names,
            fasta_file=filtered_fa,
            out_plasmid_fa=plasmid_fa,
            out_chrom_fa=chrom_fa,
            seqkit=seqkit,
        )
        logger.success("Plasmid contigs: {}", plasmid_fa)
        logger.success("Chromosome contigs: {}", chrom_fa)
    elif dry_run:
        logger.info("Dry run enabled; skip extracting and splitting contigs.")
    click.echo("finished")
                    

if __name__ == "__main__":
    main()







    
