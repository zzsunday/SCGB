############################################################################
# -*- coding:utf8 -*-
# version: 1.0
# script: FilterFalsePlasmid.py
# author: sundayzz
# mail: sundayzz@sina.cn
# This script splits assembly contigs into plasmid and chromosome FASTA files.
############################################################################

import sys
import shlex
import subprocess
from pathlib import Path
from typing import Union, Tuple
import click
from loguru import logger
from Bio import SeqIO


# =============== utils =============== #
def ensure_file(path_str: Union[str, Path]) -> str:
    """Ensure file exists (no exec bit check), return string path."""
    p = Path(path_str).expanduser()
    if p.is_file():
        return str(p)
    raise FileNotFoundError(f"File not found: {p}")


def q(value: Union[str, Path]) -> str:
    """Shell-quote a path or scalar."""
    return shlex.quote(str(value))


def run(cmd: str) -> int:
    """Run a shell command; exit the program if it fails."""
    logger.info("CMD: {}", cmd)
    try:
        subprocess.run(cmd, shell=True, check=True)
        return 0
    except subprocess.CalledProcessError as e:
        logger.error("Command failed with exit code {}: {}", e.returncode, cmd)
        sys.exit(e.returncode)


# =============== core I/O ===============
def count_fasta_records(fasta_file: Union[str, Path]) -> int:
    """Count FASTA records with Biopython."""
    with Path(fasta_file).expanduser().open("r") as fh:
        return sum(1 for _ in SeqIO.parse(fh, "fasta"))


def split_plasmid_and_chromosome(
    names_file: Union[str, Path],
    fasta_file: Union[str, Path],
    out_plasmid_fa: Union[str, Path],
    out_chrom_fa: Union[str, Path],
    seqkit: Union[str, Path] = "seqkit",
) -> Tuple[Path, Path, int, int, int]:
    """
    Split sequences with seqkit grep.

    Plasmids:
      seqkit grep -n -r -f names.txt assembly.fa -o plasmid.fa
    Chromosomes:
      seqkit grep -n -r -v -f names.txt assembly.fa -o chromosome.fa
    """
    names_path = Path(names_file).expanduser()
    fasta_path = Path(fasta_file).expanduser()
    out_plasmid_path = Path(out_plasmid_fa).expanduser()
    out_chrom_path = Path(out_chrom_fa).expanduser()
    out_plasmid_path.parent.mkdir(parents=True, exist_ok=True)
    out_chrom_path.parent.mkdir(parents=True, exist_ok=True)

    ensure_file(names_path)
    ensure_file(fasta_path)
    cmd_plasmid = (
        f"{q(seqkit)} grep -n -r -f {q(names_path)} "
        f"{q(fasta_path)} -o {q(out_plasmid_path)}"
    )
    cmd_chrom = (
        f"{q(seqkit)} grep -n -r -v -f {q(names_path)} "
        f"{q(fasta_path)} -o {q(out_chrom_path)}"
    )
    run(cmd_plasmid)
    run(cmd_chrom)

    total = count_fasta_records(fasta_path)
    p_kept = count_fasta_records(out_plasmid_path)
    c_kept = count_fasta_records(out_chrom_path)
    logger.info(
        "Seqkit splitting done: total={}, plasmid={}, chromosome={}, out_plasmid={}, out_chrom={}",
        total, p_kept, c_kept, out_plasmid_path, out_chrom_path
    )
    return out_plasmid_path, out_chrom_path, total, p_kept, c_kept


# =============== CLI ===============
@click.command(
    help="Split assembly contigs into plasmid and chromosome FASTA files.",
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option("--infile", "-i",
              required=True,
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              help="Plain text file of plasmid contig names (one per line).")
@click.option("--fa", "-a",
              required=True,
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              help="Assembly contig FASTA file.")
@click.option("--outdir", "-o",
              required=True,
              type=click.Path(file_okay=False, path_type=Path),
              help="Output directory.")
@click.option("--prefix", "-p",
              default="plasmid",
              show_default=True,
              help="Output file prefix.")
@click.option("--seqkit", default="seqkit", show_default=True,
              help="seqkit executable used to split plasmid/chromosome FASTA files.")
def main(
    infile: Path,
    fa: Path,
    outdir: Path,
    prefix: str,
    seqkit: str,
):
    """Split assembly FASTA into plasmid and chromosome FASTA files."""
    # logging
    outdir.mkdir(parents=True, exist_ok=True)
    logger.remove()
    logger.add(sys.stdout, level="INFO",
               format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level}</level> | {message}")
    logger.add(outdir / "extract.log", level="INFO",
               format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}",
               mode="a", enqueue=True)

    # split FASTA into plasmid/chromosome
    out_plasmid = outdir / f"{prefix}.plasmid_contigs.fa"
    out_chrom = outdir / f"{prefix}.chromosome_contigs.fa"
    p_fa, c_fa, total, p_kept, c_kept = split_plasmid_and_chromosome(
        names_file=infile,
        fasta_file=fa,
        out_plasmid_fa=out_plasmid,
        out_chrom_fa=out_chrom,
        seqkit=seqkit,
    )
    logger.info("Total: {} | Plasmid: {} -> {} | Chromosome: {} -> {}", total, p_kept, p_fa, c_kept, c_fa)
    click.echo(f"Plasmid FASTA: {p_fa}")
    click.echo(f"Chromosome FASTA: {c_fa}")


if __name__ == "__main__":
    main()
