# Fraction-Resolved Translational Profiling of Vincristine Resistance in DLBCL

This repository contains the R and Python analysis code used in the thesis:

**Translational Reprogramming as a Feature of Vincristine Resistance in Diffuse Large B-Cell Lymphoma**
Nimeree Muntasir

## Overview

This project uses translation-complex profiling (TCP-seq) combined with Oxford Nanopore direct RNA sequencing to characterise fraction-resolved translational regulation in vincristine-sensitive (SU-DHL8) and vincristine-resistant (SU-DHL8.R) DLBCL cell lines. TCP-seq separates ribosome-associated RNA populations into three fractions — scanning pre-initiation complexes (SSU), elongating 80S monosomes (RS), and collision-associated disomes (DS) — enabling transcriptome-wide resolution of distinct translational stages simultaneously.

## Repository Structure
config/         Path configuration for input, analysis, and output directories

data/           External reference files (see data/external/README.md)

metadata/       Sample sheet and experimental metadata

scripts/

01_metric_pipeline/       Build transcript and gene-level translation metrics from P-site data

02_differential_analysis/ Limma differential analysis for P-site fraction counts and translation metrics

03_validation/            QC and validation plots for TCP-seq fractions and translation scores

04_results/               Figure generation for all main thesis results

05_candidate_genes/       Multi-stage profiling and visualisation of representative candidate genes

archive/                  Exploratory and superseded scripts not used in final analysis

## Requirements

All primary analyses were performed in R. The following packages are required:

**R packages:**
- riboWaltz (P-site estimation)
- limma, edgeR (differential analysis)
- gprofiler2 (pathway over-representation analysis)
- ggplot2, ggrepel, patchwork, cowplot (visualisation)
- data.table, dplyr, stringr (data manipulation)
- Biostrings, GenomicRanges, rtracklayer, Rsamtools (genomic data handling)
- openxlsx (Excel output)

**Python packages:**
- pandas

Run `session_info.R` to see the exact package versions used in this analysis.

## Data Access

External reference files required by scripts in `data/external/` are described in `data/external/README.md`.

Reference genome and annotation files used:
- Homo_sapiens.GRCh38 genome assembly
- Ensembl release 114 gene annotation (Homo_sapiens.GRCh38.114.chr.gtf)
- Ensembl CDS sequences (Homo_sapiens.GRCh38.cds.all.fa.gz)

## Running the Pipeline

Scripts are numbered to reflect the intended execution order:

1. Configure input and output paths in `config/paths.R`
2. Run scripts in `01_metric_pipeline/` to build transcript and gene-level translation index matrices
3. Run scripts in `02_differential_analysis/` to perform limma analysis on P-site fraction counts and composite translation metrics
4. Run scripts in `03_validation/` to generate QC and validation figures
5. Run scripts in `04_results/` to generate all main results figures
6. Run scripts in `05_candidate_genes/` for multi-stage profiling of representative candidates

## Citation

If you use this code, please cite the associated thesis and the TCP-seq method:

Shirokikh NE, Archer SK, Beilharz TH, Powell D, Preiss T. (2017). Translation complex profile sequencing to study the in vivo dynamics of mRNA-ribosome interactions during translation initiation, elongation and termination. Nature Protocols, 12(4), 697-731.
