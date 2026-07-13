# SCGB Python CLI

SCGB is a Python command-line toolkit for single-cell-guided metagenomic binning and plasmid-host linkage.

Users provide their own metagenomic assembly FASTA. SCGB now exposes three main modules: predict assembly contigs as chromosome/plasmid contigs, bin chromosome contigs into genomes, and link plasmid contigs to binned genomes with single-cell reads.

## Installation

Create the recommended Python and tool environments:

```bash
bash install/create_envs.sh
```

If you use Conda instead of Micromamba:

```bash
bash install/create_envs.sh --manager conda
```

Activate the main SCGB environment:

```bash
micromamba activate scgb
pip install -e .
scgb --help
```

## Main Workflow

### 1. Predict assembly contigs as chromosome/plasmid contigs

```bash
scgb predict \
  --infile assembly.fa \
  --p1 genomad \
  --p2 plasmer \
  --db1 /path/to/genomad_db \
  --db2 /path/to/plasmer_db \
  --outdir 01_predict_contigs \
  --prefix sample1 \
  --thread 16 \
  --env-manager micromamba \
  --genomad-env genomad \
  --plasmer-env plsmer
```

Outputs include `sample1.plasmid_contigs.fa` and `sample1.chromosome_contigs.fa`.
Contig splitting uses `seqkit grep -n -r -f` by default; use `--seqkit /path/to/seqkit` if needed.

### 2. Bin chromosome contigs

```bash
scgb binning \
  --infile 01_predict_contigs/sample1.chromosome_contigs.fa \
  --r1 02_CleanData/mock_clean_R1.fastq \
  --r2 02_CleanData/mock_clean_R2.fastq \
  --outdir 02_bin_chromosomes \
  --prefix sample1 \
  --threads 16
```

The final binning step runs dRep dereplication by default with settings matching:

```bash
dRep dereplicate 02_bin_chromosomes/Drep --S_algorithm ANImf -nc .5 -l 10000 -N50W 0 -sizeW 1 --ignoreGenomeQuality --clusterAlg single -g "02_bin_chromosomes/binette/final_bins/*.fa" -p 32
```

Set `--drep-genomes  to use a different genome FASTA set, or `--skip-drep` to stop after Binette.


```

Inspect commands before running heavy tools:

```bash
scgb binning ... --dry-run
```

### 3. Link plasmid contigs to binned genomes with single-cell reads

```bash
scgb link \
  -fq sample1.all_single_cell.fq \
  --single-reads-dir 03_single_cell_reads \
  -c counts.csv \
  --meta barcodes.tsv \
  -fa_dir 02_bin_chromosomes/Drep/dereplicated_bacteria \
  -p 01_predict_contigs/sample1.plasmid_contigs.fa \
  --outdir 03_link_plasmids \
  --prefix sample1 \
  --threads 16
```

This module concatenates dereplicated genome FASTA files, maps the merged single-cell reads to them, computes `bc2ref.csv` and `purity.tsv`, extracts per-genome single-cell reads, maps those reads to plasmids, writes `plasmid_to_genome.tsv`, and writes genome FASTA files with linked plasmids appended under `<outdir>/genomes`.



This project is Python-CLI-only. Snakemake files are intentionally not included in this clean folder.
