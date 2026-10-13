#!/usr/bin/env python
"""Summarise the picked genomes: genome length, total genes, and genes in the COG set they were picked in.

Writes picks_summary.tsv: one row per picked genome per set.
"""
import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--picks", required=True, help="picks.tsv from pick_genomes.py")
    p.add_argument("--presence", required=True, help="Panaroo gene_presence_absence.Rtab")
    p.add_argument("--qc-report", required=True, help="qc_report.tsv from qc_filter.py (genome lengths)")
    p.add_argument("--distances", required=True, help="genome_distances.py output directory (genes per set)")
    return p.parse_args()


def main():
    args = parse_args()
    picks = pd.read_csv(args.picks, sep="\t")

    # Gene clusters present in each genome, across the whole pangenome
    presence = pd.read_csv(args.presence, sep="\t", index_col="Gene")
    total_genes = presence.sum().rename("total_genes")

    qc = pd.read_csv(args.qc_report, sep="\t", index_col="genome_id")
    length = (qc["sum_len"] / 1e6).round(2).rename("genome_length_mb")

    # Genes present in the COG set each genome was picked in
    set_genes = pd.concat([
        pd.read_csv(Path(args.distances) / f"{s}.pcoa.tsv", sep="\t", usecols=["genome_id", "genes_present"])
        .assign(set=s)
        for s in picks["set"].unique()
    ]).rename(columns={"genes_present": "set_genes"})

    summary = (picks[["set", "label", "genome_id", "phylogroup"]]
               .join(length, on="genome_id")
               .join(total_genes, on="genome_id")
               .merge(set_genes, on=["set", "genome_id"], how="left"))
    summary.to_csv("picks_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
