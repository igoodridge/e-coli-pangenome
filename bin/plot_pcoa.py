#!/usr/bin/env python
"""Plot the PCoA of each COG set: genomes coloured by phylogroup, every picked genome ringed and
labelled, next to the variance explained by each component.

Writes plots/<set>.png.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from genome_distances import load_set
from pick_genomes import load_phylogroups

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e4e3df"
BAR = "#2a78d6"
# Fixed colour per phylogroup, so a phylogroup looks the same in every plot
PHYLOGROUP_COLOURS = {
    "B2": "#2a78d6", "D": "#eb6834", "B1": "#1baf7a", "A": "#eda100", "F": "#e87ba4",
    "C": "#008300", "G": "#6250d6", "E": "#e34948", "U": "#8c5a2b", "U/cryptic": "#3f8f9f",
    "cryptic": "#9c9a3c",
}
# Phylogroups not listed above, and genomes without one
EXTRA_COLOURS = ["#7a4f9a", "#c2410c", "#4d7c0f"]
UNKNOWN = "#b9b7b0"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--distances", required=True, help="genome_distances.py output directory")
    p.add_argument("--picks", help="picks.tsv from pick_genomes.py")
    p.add_argument("--phylogroups", help="TSV with guuid and phylogroup columns")
    p.add_argument("--extra-genomes", help="extra genomes samplesheet (genome_id, fasta, phylogroup)")
    p.add_argument("--outdir", default="plots")
    return p.parse_args()


def style(ax):
    ax.set_facecolor(SURFACE)
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)
    for side in ["left", "bottom"]:
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_SECONDARY, labelsize=9)


def plot_set(coords, variance, picks=None, pg=None, title=""):
    """Figure with the PCo1/PCo2 scatter and the variance explained per component."""
    """picks: every picked genome (genome_id, label), ringed and labelled on the plot."""
    picks = pd.DataFrame(columns=["genome_id", "label"]) if picks is None else picks
    fig, (ax, ax_var) = plt.subplots(1, 2, figsize=(15, 5.5), width_ratios=[1.7, 1], facecolor=SURFACE)
    explained = variance.set_index("component")["explained"]

    # Scatter: one colour per phylogroup, largest drawn first so smaller groups stay visible on top
    group = pd.Series("unknown", index=coords.index)
    if pg is not None and not pg.empty:
        group = pg.reindex(coords.index).fillna("unknown")
    extras = iter(EXTRA_COLOURS)
    for g in group.value_counts().index:
        colour = UNKNOWN if g == "unknown" else PHYLOGROUP_COLOURS.get(g) or next(extras, UNKNOWN)
        pts = coords[group == g]
        ax.scatter(pts["PCo1"], pts["PCo2"], s=10, color=colour, alpha=0.75, linewidths=0,
                   label=f"{g} ({len(pts)})")

    labels = picks.drop_duplicates("genome_id").set_index("genome_id")["label"]
    picked = coords.loc[coords.index.intersection(labels.index)]
    ax.scatter(picked["PCo1"], picked["PCo2"], s=90, facecolors="none", edgecolors=INK, linewidths=1.2,
               label="picked", zorder=3)
    # Each label sits against its own ring, on whichever side is furthest from other points and labels
    # (in screen space, so the choice holds whatever the axis scales); a halo keeps it legible over points
    ax.autoscale_view()
    to_px = ax.transData.transform
    px_per_pt = fig.dpi / 72
    others = to_px(coords[["PCo1", "PCo2"]].to_numpy())
    corners = [(6, 3, "left", "bottom"), (6, -3, "left", "top"), (-6, 3, "right", "bottom"), (-6, -3, "right", "top")]
    placed = []
    (x0, y0), (x1, y1) = ax.get_window_extent().get_points()
    for gid, row in picked.iterrows():
        xy = to_px([(row["PCo1"], row["PCo2"])])[0]
        text = labels[gid]
        half = px_per_pt * np.array([3.5 * len(text), 6])  # rough half size of the label
        near = np.vstack([others[~np.all(np.isclose(others, xy), axis=1)], *placed])

        def centre(corner):
            dx, dy = corner[:2]
            return xy + px_per_pt * np.array([dx, dy]) + np.sign([dx, dy]) * half

        def clearance(corner):
            c = centre(corner)
            inside = x0 < c[0] - half[0] and c[0] + half[0] < x1 and y0 < c[1] - half[1] and c[1] + half[1] < y1
            return inside, np.min(np.linalg.norm(near - c, axis=1)) if len(near) else np.inf

        best = max(corners, key=clearance)
        placed.append(centre(best)[None])
        dx, dy, ha, va = best
        ax.annotate(text, (row["PCo1"], row["PCo2"]), xytext=(dx, dy), textcoords="offset points", ha=ha, va=va,
                    fontsize=10, fontweight="bold", color=INK, zorder=4,
                    path_effects=[pe.withStroke(linewidth=3, foreground=SURFACE)])

    style(ax)
    ax.set_xlabel(f"PCo1 ({explained['PCo1']:.1%})", color=INK_SECONDARY)
    ax.set_ylabel(f"PCo2 ({explained['PCo2']:.1%})", color=INK_SECONDARY)
    legend = ax.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="center left",
                       bbox_to_anchor=(1.01, 0.5), markerscale=1.8, title="Phylogroup", title_fontsize=9)
    # markerscale suits the small phylogroup dots but blows up the picked ring, so shrink that one back
    for text, handle in zip(legend.get_texts(), legend.legend_handles):
        if text.get_text() == "picked":
            handle.set_sizes([40])

    # Variance explained per component; PCo1 and PCo2 (the scatter's axes) are labelled
    ax_var.bar(explained.index, explained.values * 100, color=BAR, width=0.7)
    for i, comp in enumerate(["PCo1", "PCo2"]):
        ax_var.annotate(f"{explained[comp]:.1%}", (i, explained[comp] * 100), xytext=(0, 3),
                        textcoords="offset points", ha="center", fontsize=9, color=INK)
    style(ax_var)
    ax_var.set_ylabel("Variance explained (%)", color=INK_SECONDARY)
    ax_var.tick_params(axis="x", labelrotation=45)

    fig.suptitle(title, color=INK, fontsize=13, x=0.01, ha="left")
    fig.tight_layout()
    return fig


def main():
    args = parse_args()
    dist_dir = Path(args.distances)
    sets = pd.read_csv(dist_dir / "sets.tsv", sep="\t")
    variance = pd.read_csv(dist_dir / "variance.tsv", sep="\t")
    picks = pd.read_csv(args.picks, sep="\t") if args.picks else pd.DataFrame(columns=["set", "genome_id", "label"])
    pg = load_phylogroups(args.phylogroups, args.extra_genomes)

    out = Path(args.outdir)
    out.mkdir(exist_ok=True)
    for _, s in sets.iterrows():
        _, coords = load_set(dist_dir, s["set"])
        here = picks.loc[picks["set"] == s["set"], "label"]
        title = (f"{s['set']}: {s['genes']} genes, {s['genomes']} genomes"
                 + (f"  |  picked in this set: {', '.join(here)}" if len(here) else ""))
        fig = plot_set(coords, variance[variance["set"] == s["set"]], picks, pg, title)
        fig.savefig(out / f"{s['set']}.png", dpi=150, facecolor=SURFACE)
        plt.close(fig)
        print(f"plots/{s['set']}.png")


if __name__ == "__main__":
    main()
