from pathlib import Path
import sys
_HERE = Path(__file__).resolve()
_REPO_ROOT = next(p for p in _HERE.parents if (p / "config" / "paths.py").exists())
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
from config.paths import analysis_path, input_path, project_resource_path, external_path

import os
import sys

import pandas as pd

supp_xls = project_resource_path("HNRNPD_AUF1_targets", "Yoon2014_AUF1_PARCLIP_supp_table1.xls")
out_dir = analysis_path("Limma_translation_metrics_lfc0.7_rawP0.05", "Multi_metric_integration", "HNRNPD_AUF1_targets")
os.makedirs(out_dir, exist_ok=True)

genes = pd.read_excel(supp_xls, sheet_name="Genes", engine="xlrd")
genes.columns = [str(c).strip().replace("\n", " ") for c in genes.columns]
genes = genes.rename(
    columns={
        "TranscriptID": "transcript_id",
        "GeneName": "gene_name",
        "AUF1_target sites": "auf1_target_sites",
        "AUF1_T-to-C": "auf1_t_to_c",
        "HUR-target sites": "hur_target_sites",
        "HUR-T-to-C": "hur_t_to_c",
    }
)
genes = genes[genes["gene_name"].notna()].copy()
genes["gene_name"] = genes["gene_name"].astype(str).str.strip()
genes = genes[(genes["gene_name"] != "") & (genes["gene_name"].str.upper() != "UNKNOWN")]

for col in ["auf1_target_sites", "auf1_t_to_c", "hur_target_sites", "hur_t_to_c"]:
    genes[col] = pd.to_numeric(genes[col], errors="coerce").fillna(0)

gene_level = (
    genes.groupby("gene_name", as_index=False)
    .agg(
        n_target_transcripts=("transcript_id", "nunique"),
        total_auf1_target_sites=("auf1_target_sites", "sum"),
        total_auf1_t_to_c=("auf1_t_to_c", "sum"),
        total_hur_target_sites=("hur_target_sites", "sum"),
        total_hur_t_to_c=("hur_t_to_c", "sum"),
    )
    .sort_values(["total_auf1_target_sites", "total_auf1_t_to_c"], ascending=False)
)
gene_level["auf1_parclip_target"] = gene_level["total_auf1_target_sites"] > 0

all_path = os.path.join(out_dir, "Yoon2014_AUF1_PARCLIP_gene_level_targets.csv")
target_path = os.path.join(out_dir, "Yoon2014_AUF1_PARCLIP_gene_symbols.txt")
gene_level.to_csv(all_path, index=False)
with open(target_path, "w") as handle:
    for gene in gene_level.loc[gene_level["auf1_parclip_target"], "gene_name"]:
        handle.write(f"{gene}\n")

print("AUF1 gene-level rows:", len(gene_level))
print("AUF1 target genes:", int(gene_level["auf1_parclip_target"].sum()))
print("Top target genes:", ", ".join(gene_level.head(20)["gene_name"].tolist()))
print("Output:", out_dir)
