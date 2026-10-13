#!/usr/bin/env python
"""Jaccard distances and PCoA between genomes, for each COG gene set.

Writes into distances/:
  <set>.npz        genome x genome Jaccard distances (load with load_set)
  <set>.pcoa.tsv   PCo1-PCo10 coordinates and genes present, per genome
  variance.tsv     variance explained by each component, per set
  sets.tsv         genes and genomes in each set
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

N_COMPONENTS = 10


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--presence", required=True, help="Panaroo gene_presence_absence.Rtab")
    p.add_argument("--annotations", required=True, help="gene_annotations.tsv (Gene, COG_category)")
    p.add_argument("--cog-sets", required=True, help="COG sets YAML (see conf/cog_sets.yml)")
    p.add_argument("--outdir", default="distances")
    return p.parse_args()


def jaccard(X):
    """Pairwise Jaccard distance between the rows of a 0/1 matrix."""
    X = X.astype(np.float32)
    shared = X @ X.T
    n = X.sum(axis=1)
    union = n[:, None] + n[None, :] - shared
    return 1 - np.divide(shared, union, out=np.zeros_like(shared), where=union > 0)


def pcoa(D, n_components=N_COMPONENTS):
    """Principal coordinates of a distance matrix, and the variance each component explains."""
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D.astype(np.float64) ** 2) @ J
    vals, vecs = np.linalg.eigh(B)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    coords = vecs[:, :n_components] * np.sqrt(np.maximum(vals[:n_components], 0))
    # Conservative: negative eigenvalues count against the total
    explained = vals[:n_components] / np.abs(vals).sum()
    return coords, explained


def load_set(dist_dir, name):
    """Distance matrix (DataFrame) and PCoA table for one COG set."""
    dist_dir = Path(dist_dir)
    npz = np.load(dist_dir / f"{name}.npz")
    d = pd.DataFrame(npz["distances"], index=npz["genomes"], columns=npz["genomes"])
    coords = pd.read_csv(dist_dir / f"{name}.pcoa.tsv", sep="\t", index_col="genome_id")
    return d, coords


def main():
    args = parse_args()
    with open(args.cog_sets) as fh:
        cog_sets = yaml.safe_load(fh)

    binary = pd.read_csv(args.presence, sep="\t", index_col="Gene").T   # genomes x genes
    cats = pd.read_csv(args.annotations, sep="\t", usecols=["Gene", "COG_category"]).dropna()
    cats = cats.set_index("Gene")["COG_category"]

    out = Path(args.outdir)
    out.mkdir(exist_ok=True)
    sets, variance = [], []
    for name, letters in cog_sets.items():
        if letters == "all":
            genes = binary.columns
        else:
            genes = cats.index[cats.str.contains(f"[{''.join(letters)}]")]
        table = binary.loc[:, binary.columns.intersection(genes)]
        table = table.loc[table.sum(axis=1) > 0]   # genomes with none of these genes have no distance

        D = jaccard(table.to_numpy())
        coords, explained = pcoa(D)
        np.savez_compressed(out / f"{name}.npz", distances=D.astype(np.float32), genomes=table.index.to_numpy(str))
        pcoa_table = pd.DataFrame(coords, index=table.index.rename("genome_id"),
                                  columns=[f"PCo{i}" for i in range(1, coords.shape[1] + 1)])
        pcoa_table["genes_present"] = table.sum(axis=1)
        pcoa_table.to_csv(out / f"{name}.pcoa.tsv", sep="\t", float_format="%.6f")

        sets.append({"set": name, "genes": table.shape[1], "genomes": table.shape[0]})
        variance += [{"set": name, "component": f"PCo{i}", "explained": round(v, 5)}
                     for i, v in enumerate(explained, 1)]
        print(f"{name}: {table.shape[1]} genes, {table.shape[0]} genomes, "
              f"PCo1 {explained[0]:.1%}, PCo2 {explained[1]:.1%}")

    pd.DataFrame(sets).to_csv(out / "sets.tsv", sep="\t", index=False)
    pd.DataFrame(variance).to_csv(out / "variance.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
