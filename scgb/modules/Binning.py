#!/usr/bin/env python3
############################################################################
# -*- coding:utf-8 -*-
# script: Binning.py
# Wrapper for MAG binning with COMEBin, SemiBin2, MetaDecoder, Binette and dRep.
############################################################################

from __future__ import annotations

import gzip
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Union

import click
from loguru import logger


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


def shell_cmd(env_manager: str, env_name: str | None) -> str:
    """Return bash, optionally wrapped by an environment manager."""
    return tool_cmd("bash", env_manager, env_name)


def ensure_file(path_str: Union[str, Path]) -> Path:
    """Ensure a file exists and return it as a Path."""
    path = Path(path_str).expanduser()
    if not path.is_file():
        raise FileNotFoundError(f"File not found: {path}")
    return path


def ensure_executable(path_str: Union[str, Path]) -> Path:
    """Ensure an executable exists and return it as a Path."""
    path = ensure_file(path_str)
    if not path.exists():
        raise FileNotFoundError(f"Executable not found: {path}")
    return path


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


def reset_dir(path: Path, overwrite: bool) -> None:
    """Create an output directory, optionally replacing an old one."""
    if path.exists() and overwrite:
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def gunzip_outputs(directory: Path, dry_run: bool = False) -> None:
    """Decompress all .gz files under a directory and remove the compressed copies."""
    if dry_run:
        logger.info("Dry run enabled; would gunzip files under {}", directory)
        return
    if not directory.exists():
        logger.warning("Skip gunzip; directory does not exist: {}", directory)
        return

    gz_files = sorted(directory.rglob("*.gz"))
    if not gz_files:
        logger.info("No .gz files found under {}", directory)
        return

    for gz_file in gz_files:
        out_file = gz_file.with_suffix("")
        logger.info("Gunzip: {} -> {}", gz_file, out_file)
        with gzip.open(gz_file, "rb") as fin, out_file.open("wb") as fout:
            shutil.copyfileobj(fin, fout)
        gz_file.unlink()


def build_bowtie2_bam(
    assembly: Path,
    r1: Path,
    r2: Path,
    outdir: Path,
    prefix: str,
    threads: int,
    bowtie2_build: Union[str, Path],
    bowtie2: Union[str, Path],
    samtools: Union[str, Path],
    dry_run: bool = False,
) -> Path:
    """
    Build a Bowtie2 index, align paired reads, and produce a sorted indexed BAM.

    This BAM is reused by both SemiBin2 and MetaDecoder.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    index_prefix = outdir / assembly.name
    sam = outdir / f"{prefix}.sam"
    bam = outdir / f"{prefix}.bam"
    mapped_bam = outdir / f"{prefix}.mapped.bam"
    sorted_bam = outdir / f"{prefix}.mapped.sorted.bam"

    run(f"{bowtie2_build} -f {q(assembly)} {q(index_prefix)}", dry_run=dry_run)
    run(
        f"{bowtie2} -q --fr -x {q(index_prefix)} -1 {q(r1)} -2 {q(r2)} "
        f"-S {q(sam)} -p {threads}",
        dry_run=dry_run,
    )
    run(f"{samtools} view -h -b -S {q(sam)} -o {q(bam)}", dry_run=dry_run)
    run(f"{samtools} view -b -F 4 {q(bam)} -o {q(mapped_bam)}", dry_run=dry_run)
    run(f"{samtools} sort {q(mapped_bam)} -o {q(sorted_bam)}", dry_run=dry_run)
    run(f"{samtools} index {q(sorted_bam)}", dry_run=dry_run)
    return sorted_bam


def run_semibin2(
    semibin2: Union[str, Path],
    assembly: Path,
    sorted_bam: Path,
    outdir: Path,
    threads: int,
    dry_run: bool = False,
) -> Path:
    """Run SemiBin2 single_easy_bin."""
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = (
        f"{semibin2} single_easy_bin "
        f"-i {q(assembly)} -b {q(sorted_bam)} -o {q(outdir)} --threads {threads}"
    )
    run(cmd, dry_run=dry_run)
    bin_dir = outdir / "output_bins"
    gunzip_outputs(bin_dir, dry_run=dry_run)
    return bin_dir


def run_metadecoder(
    metadecoder: Union[str, Path],
    assembly: Path,
    sorted_bam: Path,
    outdir: Path,
    prefix: str,
    threads: int,
    dry_run: bool = False,
) -> Path:
    """Run MetaDecoder coverage, seed and cluster."""
    outdir.mkdir(parents=True, exist_ok=True)
    coverage = outdir / f"{prefix}.coverage"
    seed = outdir / f"{prefix}.seed"
    cluster_prefix = outdir / prefix

    run(
        f"{metadecoder} coverage --threads {threads} -b {q(sorted_bam)} -o {q(coverage)}",
        dry_run=dry_run,
    )
    run(
        f"{metadecoder} seed --threads {threads} -f {q(assembly)} -o {q(seed)}",
        dry_run=dry_run,
    )
    run(
        f"{metadecoder} cluster -f {q(assembly)} -c {q(coverage)} "
        f"-s {q(seed)} -o {q(cluster_prefix)}",
        dry_run=dry_run,
    )
    return outdir


def run_comebin(
    assembly: Path,
    r1: Path,
    r2: Path,
    outdir: Path,
    prefix: str,
    threads: int,
    env_manager: str,
    comebin_env: str,
    forward_suffix: str,
    reverse_suffix: str,
    dry_run: bool = False,
) -> Path:
    """
    Run COMEBin coverage generation and COMEBin binning.

    COMEBin helper scripts are found inside the COMEBin environment.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    cov_out = outdir
    cov_work_files = cov_out / "work_files"
    bin_out = outdir / f"{prefix}_comebin"
    cov_out.mkdir(parents=True, exist_ok=True)

    bash = shell_cmd(env_manager, comebin_env)
    cov_script = (
        'GEN_COV="$(find "$CONDA_PREFIX" -name gen_cov_file.sh | head -n 1)"; '
        'if [ -z "$GEN_COV" ]; then echo "gen_cov_file.sh not found in $CONDA_PREFIX" >&2; exit 1; fi; '
        f'bash "$GEN_COV" -a {q(assembly)} -t {threads} '
        f'-f {q(forward_suffix)} -r {q(reverse_suffix)} -o {q(cov_out)} {q(r1)} {q(r2)}'
    )
    bin_script = (
        'RUN_COMEBIN="$(find "$CONDA_PREFIX" -name run_comebin.sh | head -n 1)"; '
        'if [ -z "$RUN_COMEBIN" ]; then echo "run_comebin.sh not found in $CONDA_PREFIX" >&2; exit 1; fi; '
        f'if [ ! -d {q(cov_work_files)} ]; then echo "COMEBin work_files not found: {q(cov_work_files)}" >&2; exit 1; fi; '
        f'bash "$RUN_COMEBIN" -a {q(assembly)} -p {q(cov_work_files)} '
        f'-o {q(bin_out)} -t {threads}'
    )

    run(
        f"{bash} -lc {q(cov_script)}",
        dry_run=dry_run,
    )
    run(
        f"{bash} -lc {q(bin_script)}",
        dry_run=dry_run,
    )
    return bin_out / "comebin_res" / "comebin_res_bins"


def run_binette(
    binette: Union[str, Path],
    assembly: Path,
    semibin_bins: Path,
    metadecoder_bins: Path,
    comebin_bins: Path,
    outdir: Path,
    threads: int,
    dry_run: bool = False,
) -> Path:
    """Refine bins from SemiBin2, MetaDecoder and COMEBin with Binette."""
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = (
        f"{binette} --bin_dirs {q(semibin_bins)} {q(metadecoder_bins)} {q(comebin_bins)} "
        f"--contigs {q(assembly)} --threads {threads} --outdir {q(outdir)}"
    )
    run(cmd, dry_run=dry_run)
    return outdir


def run_drep(
    drep: Union[str, Path],
    genomes: str,
    outdir: Path,
    threads: int,
    env_manager: str,
    drep_env: str,
    s_algorithm: str = "ANImf",
    nc: float = 0.5,
    min_length: int = 10000,
    n50_weight: int = 0,
    size_weight: int = 1,
    cluster_alg: str = "single",
    ignore_genome_quality: bool = True,
    extra_args: str = "",
    dry_run: bool = False,
) -> Path:
    """Dereplicate refined bins with dRep."""
    outdir.mkdir(parents=True, exist_ok=True)
    drep_cmd = tool_cmd(drep, env_manager, drep_env)
    cmd = (
        f"{drep_cmd} dereplicate {q(outdir)} "
        f"--S_algorithm {q(s_algorithm)} -nc {nc} -l {min_length} "
        f"-N50W {n50_weight} -sizeW {size_weight} "
        f"--clusterAlg {q(cluster_alg)} "
        f"-g {genomes} -p {threads}"
    )
    if ignore_genome_quality:
        cmd += " --ignoreGenomeQuality"
    if extra_args:
        cmd += f" {extra_args}"
    run(cmd, dry_run=dry_run)
    return outdir


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option("--infile", "-i", required=True,
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              help="Chromosome contig FASTA used for MAG binning.")
@click.option("--r1", "-1", required=True,
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              help="Clean paired-end read 1 FASTQ.")
@click.option("--r2", "-2", required=True,
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              help="Clean paired-end read 2 FASTQ.")
@click.option("--outdir", "-o", required=True,
              type=click.Path(file_okay=False, path_type=Path),
              help="Output directory.")
@click.option("--prefix", "-p", required=True, help="Output prefix.")
@click.option("--threads", "-t", default=16, show_default=True,
              type=click.IntRange(min=1), help="CPU threads.")
@click.option("--forward-suffix", "comebin_forward_suffix", default="_1.fastq", show_default=True,
              help="Forward read suffix passed to COMEBin gen_cov_file.sh -f.")
@click.option("--reverse-suffix", "comebin_reverse_suffix", default="_2.fastq", show_default=True,
              help="Reverse read suffix passed to COMEBin gen_cov_file.sh -r.")
@click.option("--metadecoder", default="metadecoder", show_default=True,
              help="MetaDecoder executable.")
@click.option("--semibin2", default="SemiBin2", show_default=True,
              help="SemiBin2 executable.")
@click.option("--binette", default="binette", show_default=True,
              help="Binette executable.")
@click.option("--drep", default="dRep", show_default=True,
              help="dRep executable.")
@click.option("--samtools", default="samtools", show_default=True,
              help="samtools executable.")
@click.option("--env-manager", default="conda", show_default=True,
              help="Environment manager used to run external tools.")
@click.option("--semibin-env", default="scgb_semibin", show_default=True,
              help="Environment containing SemiBin2, Bowtie2 and samtools. Empty string = use active PATH.")
@click.option("--metadecoder-env", default="scgb_metadecoder", show_default=True,
              help="Environment containing MetaDecoder. Empty string = use active PATH.")
@click.option("--comebin-env", default="scgb_comebin", show_default=True,
              help="Environment containing COMEBin. Empty string = use active PATH.")
@click.option("--binette-env", default="scgb_binette", show_default=True,
              help="Environment containing Binette. Empty string = use active PATH.")
@click.option("--drep-env", default="scgb_drep", show_default=True,
              help="Environment containing dRep. Empty string = use active PATH.")
@click.option("--drep-genomes", default=None,
              help="Genome FASTA glob for dRep -g. Default: <outdir>/binette/final_bins/*.fa")
@click.option("--drep-outdir", default="Drep", show_default=True,
              help="dRep output subdirectory inside --outdir, or an absolute path.")
@click.option("--drep-s-algorithm", default="ANImf", show_default=True,
              help="dRep secondary clustering algorithm.")
@click.option("--drep-nc", default=0.5, show_default=True, type=float,
              help="dRep minimum alignment coverage passed to -nc.")
@click.option("--drep-min-length", default=10000, show_default=True, type=click.IntRange(min=0),
              help="dRep minimum genome length passed to -l.")
@click.option("--drep-n50-weight", default=0, show_default=True, type=int,
              help="dRep N50 weight passed to -N50W.")
@click.option("--drep-size-weight", default=1, show_default=True, type=int,
              help="dRep size weight passed to -sizeW.")
@click.option("--drep-cluster-alg", default="single", show_default=True,
              help="dRep clustering algorithm passed to --clusterAlg.")
@click.option("--drep-extra-args", default="", show_default=True,
              help="Additional raw arguments appended to dRep.")
@click.option("--skip-drep", is_flag=True,
              help="Skip dRep dereplication after Binette.")
@click.option("--dry-run", is_flag=True, help="Print commands without executing them.")
@click.option("--overwrite", is_flag=True, help="Replace an existing output directory.")
def main(
    infile: Path,
    r1: Path,
    r2: Path,
    outdir: Path,
    prefix: str,
    threads: int,
    comebin_forward_suffix: str,
    comebin_reverse_suffix: str,
    metadecoder: str,
    semibin2: str,
    binette: str,
    drep: str,
    samtools: str,
    env_manager: str,
    semibin_env: str,
    metadecoder_env: str,
    comebin_env: str,
    binette_env: str,
    drep_env: str,
    drep_genomes: str | None,
    drep_outdir: str,
    drep_s_algorithm: str,
    drep_nc: float,
    drep_min_length: int,
    drep_n50_weight: int,
    drep_size_weight: int,
    drep_cluster_alg: str,
    drep_extra_args: str,
    skip_drep: bool,
    dry_run: bool,
    overwrite: bool,
) -> None:
    """
    Bin chromosome contigs with COMEBin, SemiBin2, MetaDecoder, Binette and dRep.
    """
    reset_dir(outdir, overwrite=overwrite)
    logger.remove()
    logger.add(sys.stdout, level="INFO",
               format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level}</level> | {message}")
    logger.add(outdir / "binning.log", level="INFO",
               format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}", mode="a", enqueue=True)

    assembly = ensure_file(infile)
    read1 = ensure_file(r1)
    read2 = ensure_file(r2)

    logger.info("Assembly: {}", assembly)
    logger.info("Reads: {} | {}", read1, read2)
    logger.info("Output directory: {}", outdir)
    logger.info("Environment manager: {}", env_manager)
    logger.info("Tool environments: semibin={} metadecoder={} comebin={} binette={} drep={}",
                semibin_env or "PATH", metadecoder_env or "PATH",
                comebin_env or "PATH", binette_env or "PATH", drep_env or "PATH")

    sembin_dir = outdir / "sembin"
    metadecoder_dir = outdir / "metadecoder"
    comebin_dir = outdir / "comebin"
    binette_dir = outdir / "binette"
    drep_dir = Path(drep_outdir).expanduser()
    if not drep_dir.is_absolute():
        drep_dir = outdir / drep_dir

    bowtie2_build_cmd = tool_cmd("bowtie2-build", env_manager, semibin_env)
    bowtie2_cmd = tool_cmd("bowtie2", env_manager, semibin_env)
    samtools_cmd = tool_cmd(samtools, env_manager, semibin_env)
    semibin2_cmd = tool_cmd(semibin2, env_manager, semibin_env)
    metadecoder_cmd = tool_cmd(metadecoder, env_manager, metadecoder_env)
    binette_cmd = tool_cmd(binette, env_manager, binette_env)

    logger.info("===== Step 1: Build Bowtie2 BAM for SemiBin2 and MetaDecoder =====")
    sorted_bam = build_bowtie2_bam(
        assembly=assembly,
        r1=read1,
        r2=read2,
        outdir=sembin_dir,
        prefix=f"{prefix}_sembin",
        threads=threads,
        bowtie2_build=bowtie2_build_cmd,
        bowtie2=bowtie2_cmd,
        samtools=samtools_cmd,
        dry_run=dry_run,
    )

    logger.info("===== Step 2: Run SemiBin2 =====")
    semibin_bins = run_semibin2(
        semibin2=semibin2_cmd,
        assembly=assembly,
        sorted_bam=sorted_bam,
        outdir=sembin_dir / "output",
        threads=threads,
        dry_run=dry_run,
    )

    logger.info("===== Step 3: Run MetaDecoder =====")
    metadecoder_bins = run_metadecoder(
        metadecoder=metadecoder_cmd,
        assembly=assembly,
        sorted_bam=sorted_bam,
        outdir=metadecoder_dir,
        prefix=f"{prefix}_metadecoder",
        threads=threads,
        dry_run=dry_run,
    )

    logger.info("===== Step 4: Run COMEBin =====")
    comebin_bins = run_comebin(
        assembly=assembly,
        r1=read1,
        r2=read2,
        outdir=comebin_dir,
        prefix=prefix,
        threads=threads,
        env_manager=env_manager,
        comebin_env=comebin_env,
        forward_suffix=comebin_forward_suffix,
        reverse_suffix=comebin_reverse_suffix,
        dry_run=dry_run,
    )

    logger.info("===== Step 5: Refine bins with Binette =====")
    run_binette(
        binette=binette_cmd,
        assembly=assembly,
        semibin_bins=semibin_bins,
        metadecoder_bins=metadecoder_bins,
        comebin_bins=comebin_bins,
        outdir=binette_dir,
        threads=threads,
        dry_run=dry_run,
    )

    drep_output: Path | None = None
    if skip_drep:
        logger.info("Skip dRep dereplication because --skip-drep was set.")
    else:
        logger.info("===== Step 6: Dereplicate bins with dRep =====")
        final_bins_dir = Path(binette_dir) / "final_bins"
        genomes = drep_genomes or f"{q(final_bins_dir)}/*.fa"
        drep_output = run_drep(
            drep=drep,
            genomes=genomes,
            outdir=drep_dir,
            threads=threads,
            env_manager=env_manager,
            drep_env=drep_env,
            s_algorithm=drep_s_algorithm,
            nc=drep_nc,
            min_length=drep_min_length,
            n50_weight=drep_n50_weight,
            size_weight=drep_size_weight,
            cluster_alg=drep_cluster_alg,
            ignore_genome_quality=True,
            extra_args=drep_extra_args,
            dry_run=dry_run,
        )

    logger.success("Binning workflow finished.")
    logger.info("SemiBin2 bins: {}", semibin_bins)
    logger.info("MetaDecoder bins: {}", metadecoder_bins)
    logger.info("COMEBin bins: {}", comebin_bins)
    logger.info("Binette output: {}", binette_dir)
    if drep_output:
        logger.info("dRep output: {}", drep_output)


if __name__ == "__main__":
    main()
