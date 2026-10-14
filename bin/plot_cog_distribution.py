#!/usr/bin/env python
"""Plot the number of gene clusters in each COG category, split into core (>= 95% of genomes,
Panaroo core + soft core) and accessory, from the pipeline's gene_annotations.tsv.

A cluster with several category letters counts once in each of its categories. The poorly
characterised categories (R, S) are left out. Categories are sorted by total clusters, largest first.
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

SURFACE = "#ffffff"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e4e3df"
# Dark and light shades of one hue, so core and accessory read as parts of one total
# (and don't reuse the PCoA phylogroup colours)
CORE = "#1d4e5f"
ACCESSORY = "#8cc5c9"

# COG 2020 categories in the order of the methods table, grouped by broad function
GROUPS = {
    "Information storage and processing": {
        "J": "Translation, ribosomal structure and biogenesis", "A": "RNA processing and modification",
        "K": "Transcription", "L": "Replication, recombination and repair",
        "B": "Chromatin structure and dynamics"},
    "Cellular processes and signalling": {
        "D": "Cell cycle control, cell division, chromosome partitioning", "Y": "Nuclear structure",
        "V": "Defence mechanisms", "T": "Signal transduction mechanisms",
        "M": "Cell wall/membrane/envelope biogenesis", "N": "Cell motility", "Z": "Cytoskeleton",
        "W": "Extracellular structures", "U": "Intracellular trafficking, secretion and vesicular transport",
        "O": "Posttranslational modification, protein turnover, chaperones",
        "X": "Mobilome: prophages, transposons"},
    "Metabolism": {
        "C": "Energy production and conversion", "G": "Carbohydrate transport and metabolism",
        "E": "Amino acid transport and metabolism", "F": "Nucleotide transport and metabolism",
        "H": "Coenzyme transport and metabolism", "I": "Lipid transport and metabolism",
        "P": "Inorganic ion transport and metabolism",
        "Q": "Secondary metabolites biosynthesis, transport and catabolism"},
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--annotations", required=True, help="gene_annotations.tsv from SUMMARISE_ANNOTATIONS")
    p.add_argument("--out", required=True, help="output image (.png or .pdf)")
    return p.parse_args()


def main():
    args = parse_args()
    genes = pd.read_csv(args.annotations, sep="\t", usecols=["COG_category", "pangenome_category"]).dropna()
    core = genes["pangenome_category"].isin(["core", "soft_core"])

    rows = []
    for group, cats in GROUPS.items():
        for letter, name in cats.items():
            has = genes["COG_category"].str.contains(letter)
            if has.any():   # B and Y (eukaryotic) have no clusters
                rows.append({"group": group, "label": f"{letter}  {name}",
                             "core": (has & core).sum(), "accessory": (has & ~core).sum()})
    counts = pd.DataFrame(rows)
    counts["total"] = counts["core"] + counts["accessory"]
    counts = counts.sort_values("total", ascending=False).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8, 7.5), facecolor=SURFACE)
    y = range(len(counts))[::-1]   # first category at the top
    ax.barh(y, counts["core"], color=CORE, height=0.7,
            edgecolor=SURFACE, linewidth=1)
    ax.barh(y, counts["accessory"], left=counts["core"], color=ACCESSORY, height=0.7,
            edgecolor=SURFACE, linewidth=1)
    for yi, total in zip(y, counts["total"]):
        ax.text(total + 60, yi, f"{total:,}", va="center", fontsize=8, color=INK_SECONDARY)

    ax.set_yticks(list(y), counts["label"], fontsize=8, color=INK)

    ax.set_xlabel("Gene clusters", fontsize=9, color=INK_SECONDARY)
    ax.tick_params(axis="x", labelsize=8, colors=INK_SECONDARY)
    ax.tick_params(axis="y", length=0)
    for side in ["top", "right", "left"]:
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.set_xlim(0, counts["total"].max() * 1.12)
    ax.margins(y=0.01)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (CORE, ACCESSORY)]
    ax.legend(handles, ["Core (≥ 95% of genomes)", "Accessory (< 95% of genomes)"],
              frameon=False, fontsize=8, labelcolor=INK_SECONDARY, loc="lower right")

    fig.tight_layout()
    fig.savefig(args.out, dpi=200, facecolor=SURFACE)


if __name__ == "__main__":
    main()
