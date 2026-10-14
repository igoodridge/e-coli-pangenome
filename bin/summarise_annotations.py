#!/usr/bin/env python
"""Combine the Panaroo pangenome with eggNOG-mapper annotations and COG definitions.

Writes gene_annotations.tsv: one row per gene cluster with its frequency, pangenome category and
eggNOG + COG annotation. Its Gene column matches the rows of Panaroo's gene_presence_absence.Rtab.
"""
import argparse

import pandas as pd

# Panaroo's metadata columns in gene_presence_absence.csv; every other column is a genome
PANAROO_META = ["Gene", "Non-unique Gene name", "Annotation"]

EGGNOG_COLS = [
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "tax_ceiling",
    "farthest_donor_lineage", "COG_category", "Preferred_name",
    "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
    "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy",
    "BiGG_Reaction", "PFAMs", "annotation_confidence",
]

COG_COLS = ["COG_ID", "COG_category", "COG_description", "COG_gene", "COG_pathway", "pubmed_id", "pdb_id"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--presence", required=True, help="Panaroo gene_presence_absence.csv")
    p.add_argument("--eggnog", required=True, help="eggNOG-mapper .emapper.annotations")
    p.add_argument("--cog-defs", required=True, help="NCBI COG definitions (cog-20.def.tab)")
    return p.parse_args()


def pangenome_category(freq):
    # Same thresholds as Panaroo's summary_statistics.txt
    if freq >= 0.99:
        return "core"
    if freq >= 0.95:
        return "soft_core"
    if freq >= 0.15:
        return "shell"
    return "cloud"


def main():
    args = parse_args()

    presence = pd.read_csv(args.presence, dtype=str)
    genomes = [c for c in presence.columns if c not in PANAROO_META]

    eggnog = pd.read_csv(args.eggnog, sep="\t", comment="#", header=None, names=EGGNOG_COLS,
                         dtype=str, encoding="latin-1", index_col=False)
    # Panaroo's pan_genome_reference.fa is named by cluster, so the eggNOG query is the Panaroo Gene
    eggnog = eggnog.rename(columns={"query": "Gene", "COG_category": "COG_ID"}).replace("-", pd.NA)

    cog = pd.read_csv(args.cog_defs, sep="\t", header=None, names=COG_COLS, dtype=str, encoding="latin-1")

    # eggNOG's COG column holds a COG ID (COG0583) when there is one, otherwise category letters (S)
    annotation = eggnog.merge(cog[["COG_ID", "COG_category", "COG_description"]], on="COG_ID", how="left")
    is_letters = annotation["COG_ID"].notna() & ~annotation["COG_ID"].str.startswith("COG", na=False)
    annotation.loc[is_letters, "COG_category"] = annotation.loc[is_letters, "COG_ID"]
    annotation.loc[is_letters, "COG_ID"] = pd.NA

    present = presence[genomes].notna()
    genes = presence[PANAROO_META].assign(
        n_genomes=present.sum(axis=1),
        frequency=present.mean(axis=1),
    )
    # Categorise before rounding: e.g. 428/2854 = 0.14996 is cloud, but rounds to 0.15 (shell)
    genes["pangenome_category"] = genes["frequency"].map(pangenome_category)
    genes["frequency"] = genes["frequency"].round(4)
    genes = genes.merge(annotation, on="Gene", how="left", validate="one_to_one")
    genes.to_csv("gene_annotations.tsv", sep="\t", index=False)

    annotated = genes["seed_ortholog"].notna()
    print(f"{len(genes)} gene clusters, {len(genomes)} genomes")
    print(f"{annotated.mean():.1%} of clusters have an eggNOG hit, "
          f"{genes['COG_category'].notna().mean():.1%} have a COG category")
    print(genes["pangenome_category"].value_counts().to_string())


if __name__ == "__main__":
    main()
