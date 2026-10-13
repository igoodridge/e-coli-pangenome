#!/usr/bin/env python
"""Pick the most distant genomes in each COG set, from genome_distances.py output.

In each set, picks the 3 genomes whose closest pair is as far apart as possible (exact).
With a fixed genome, picks the 3 others that maximise the smallest distance among all 4.

Writes picks.tsv: one row per picked genome per set, labelled G1, G2... by how many sets picked it
(the fixed genome and extra genomes keep their own ID, e.g. K-12).
"""
import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from genome_distances import load_set


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--distances", required=True, help="genome_distances.py output directory")
    p.add_argument("--config", required=True, help="picking YAML (see conf/picking.yml)")
    p.add_argument("--phylogroups", help="TSV with guuid and phylogroup columns")
    p.add_argument("--extra-genomes", help="extra genomes samplesheet (genome_id, fasta, phylogroup)")
    return p.parse_args()


def pick3(D, groups=None):
    """Indices of the trio maximising its smallest pairwise distance, and that distance.

    groups: optional list of three index arrays, one pick from each.
    Pairs are tried from most to least distant, stopping once no pair can beat the best trio.
    """
    if groups is None:
        a = b = c = np.arange(len(D))
    else:
        a, b, c = groups
    best, trio = -1.0, None
    ab = D[np.ix_(a, b)]
    for flat in np.argsort(ab, axis=None)[::-1]:
        i, j = np.unravel_index(flat, ab.shape)
        if ab[i, j] <= best:
            break
        x, y = a[i], b[j]
        if x == y:
            continue
        m = np.minimum(D[x, c], D[y, c])
        m[(c == x) | (c == y)] = -1
        k = m.argmax()
        val = min(ab[i, j], m[k])
        if val > best:
            best, trio = val, (x, y, c[k])
    return (list(trio), best) if trio else ([], np.nan)


def pick(d, cfg, pg):
    """Picks for one set's distance matrix d, after min_genes/exclude filtering. Returns (genome IDs, score)."""
    fixed = cfg.get("fixed")
    candidates = d.index.difference([fixed]) if fixed else d.index
    D = d.loc[candidates, candidates].to_numpy()

    if fixed:
        # A trio's smallest distance on this matrix equals the smallest distance among the trio + fixed genome
        to_fixed = d.loc[fixed, candidates].to_numpy()
        D = np.minimum(D, np.minimum.outer(to_fixed, to_fixed))

    groups = None
    if cfg.get("phylogroups"):
        labels = pg.reindex(candidates).to_numpy()
        groups = [np.flatnonzero(labels == g) for g in cfg["phylogroups"]]
        missing = [g for g, ix in zip(cfg["phylogroups"], groups) if len(ix) == 0]
        if missing:
            raise ValueError(f"no candidate genomes in phylogroup(s) {missing}")

    trio, score = pick3(D, groups)
    picked = list(candidates[trio])
    return ([fixed] + picked if fixed else picked), score


def pick_set(d, genes_present, cfg, pg=None):
    """Picks for one COG set, applying min_genes and exclude. Returns (genome IDs, min distance, candidates)."""
    pg = pd.Series(dtype=str) if pg is None else pg
    keep = genes_present.index[genes_present >= cfg.get("min_genes", 0)].difference(cfg.get("exclude") or [])
    fixed = cfg.get("fixed")
    if fixed and fixed not in keep:
        raise ValueError(f"fixed genome {fixed} has too few genes or is excluded")
    if len(keep) < (4 if fixed else 3):
        raise ValueError(f"only {len(keep)} genomes with >= {cfg.get('min_genes', 0)} genes")
    genomes, score = pick(d.loc[keep, keep], cfg, pg)
    return genomes, score, len(keep)


def load_phylogroups(path=None, extra_genomes=None):
    """Phylogroup per genome ID, from the phylogroups TSV and the extra genomes samplesheet."""
    pg = pd.Series(dtype=str)
    if path:
        pg = pd.read_csv(path, sep="\t").drop_duplicates("guuid").set_index("guuid")["phylogroup"]
    if extra_genomes:
        extra = pd.read_csv(extra_genomes, sep="\t").set_index("genome_id")["phylogroup"].dropna()
        pg = pd.concat([pg.drop(extra.index, errors="ignore"), extra])
    return pg


def main():
    args = parse_args()
    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    pg = load_phylogroups(args.phylogroups, args.extra_genomes)
    if cfg.get("phylogroups") and pg.empty:
        raise SystemExit("phylogroups set in the config but no --phylogroups file given")

    dist_dir = Path(args.distances)
    available = pd.read_csv(dist_dir / "sets.tsv", sep="\t")["set"].tolist()
    sets = cfg.get("sets") or available
    unknown = set(sets) - set(available)
    if unknown:
        raise SystemExit(f"unknown COG set(s) {sorted(unknown)}; available: {available}")

    rows = []
    for name in sets:
        d, coords = load_set(dist_dir, name)
        try:
            genomes, score, candidates = pick_set(d, coords["genes_present"], cfg, pg)
        except ValueError as e:
            print(f"{name}: {e}, skipping")
            continue
        for g in genomes:
            rows.append({"set": name, "genome_id": g, "fixed": g == cfg.get("fixed"),
                         "min_distance": round(float(score), 5), "candidates": candidates,
                         "phylogroup": pg.get(g, pd.NA)})

    picks = pd.DataFrame(rows)
    if picks.empty:
        raise SystemExit("no picks made: check sets, min_genes, exclude and phylogroups")

    # The fixed genome and extra genomes (e.g. K-12) keep their own ID as a label;
    # the rest are labelled G1, G2... with the most frequently picked first
    named = set(pd.read_csv(args.extra_genomes, sep="\t")["genome_id"]) if args.extra_genomes else set()
    named |= {cfg["fixed"]} if cfg.get("fixed") else set()
    counts = Counter(g for g in picks["genome_id"] if g not in named)
    order = sorted(counts, key=lambda g: (-counts[g], g))
    labels = {g: f"G{i}" for i, g in enumerate(order, 1)} | {g: g for g in named}
    picks.insert(1, "label", picks["genome_id"].map(labels))
    picks.to_csv("picks.tsv", sep="\t", index=False)

    for name, grp in picks.groupby("set", sort=False):
        print(f"{name}: {', '.join(grp['label'])}  (min distance {grp['min_distance'].iloc[0]:.3f})")


if __name__ == "__main__":
    main()
