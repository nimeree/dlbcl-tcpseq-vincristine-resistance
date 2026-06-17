# External Reference Files

Scripts in this repository depend on the following external data files which
are not included in the repository due to size or licensing restrictions.

## Han et al. 2020 — Disome collision reference

**File:** NIHMS1591487-supplement-Table_S1.csv
**Source:** Han et al. (2020). Genome-wide survey of ribosome collision. Cell Reports, 31(5), 107610.
**URL:** https://doi.org/10.1016/j.celrep.2020.107610
**Used by:** 09_plot_han_2020_collision_validation.R

Download Supplementary Table S1 from the above DOI and place it in this directory.

## Yoon et al. 2014 — AUF1 PAR-CLIP targets

**File:** Yoon2014_AUF1_PARCLIP_supp_table1.xls
**Source:** Yoon et al. (2014). Functional interactions among microRNAs and long noncoding RNAs. Science, 346(6209).
**URL:** https://doi.org/10.1126/science.1257493
**Used by:** archive/29_hnrnpd_auf1_target_interaction_test.py, archive/30_test_hnrnpd_auf1_targets_in_limma_interaction.R

Download Supplementary Table 1 from the above DOI and place it in this directory.

## Ensembl Reference Files

The following Ensembl GRCh38 release 114 files are required and should be
placed in the input directory defined in config/paths.R:

- Homo_sapiens.GRCh38.114.chr.gtf
- Homo_sapiens.GRCh38.cds.all.fa.gz

Download from: https://ftp.ensembl.org/pub/release-114/
