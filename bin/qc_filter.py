#!/usr/bin/env python
"""Select genomes passing QC from seqkit stats, CheckM2 and skani reports.

Writes:
  qc_report.tsv  every genome, its metrics, qc_pass and the reasons it failed
  qc_passed.tsv  genomes passing every rule (input list for the pipeline)
"""
import argparse
import re
from pathlib import Path

import pandas as pd
import yaml


def genome_id(path):
    # Same rule as the pipeline: file name without .fasta/.fna/.fa(.gz), so IDs match across reports
    return re.sub(r"\.(fasta|fna|fa)(\.gz)?$", "", Path(path).name)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seqkit", required=True, help="seqkit stats -a -T output")
    p.add_argument("--checkm2", required=True, help="CheckM2 quality_report.tsv")
    p.add_argument("--skani", required=True, help="skani search output vs species database")

    p.add_argument("--config", required=True, help="QC thresholds YAML (see conf/qc.yml)")
    return p.parse_args()


def skani_top_hits(path, min_af):
    skani = pd.read_csv(path, sep="\t").assign(genome_id=lambda d: d["Query_file"].map(genome_id))
    hits = skani[(skani["Align_fraction_query"] >= min_af) & (skani["Align_fraction_ref"] >= min_af)]
    return (
        hits.sort_values("ANI", ascending=False)
        .drop_duplicates(subset="genome_id", keep="first")
        .set_index("genome_id")
        [["Ref_name", "ANI", "Align_fraction_query", "Align_fraction_ref"]]
        .rename(columns={"Ref_name": "top_hit_species", "ANI": "top_hit_ani",
                         "Align_fraction_query": "top_hit_af_query", "Align_fraction_ref": "top_hit_af_ref"})
    )


def main():
    args = parse_args()
    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    seqkit = (pd.read_csv(args.seqkit, sep="\t")
              .assign(genome_id=lambda d: d["file"].map(genome_id))
              .set_index("genome_id")[["file", "num_seqs", "sum_len", "N50", "GC(%)"]])
    checkm2 = (pd.read_csv(args.checkm2, sep="\t")
               .assign(genome_id=lambda d: d["Name"].map(genome_id))
               .set_index("genome_id")[["Completeness", "Contamination"]])
    skani = skani_top_hits(args.skani, cfg["skani"]["min_af"])

    ids = seqkit.index.union(checkm2.index)
    report = (pd.DataFrame(index=ids.rename("genome_id"))
              .join(seqkit).join(checkm2).join(skani))

    species = report["top_hit_species"].fillna("").str.contains(cfg["skani"]["species"], case=False, regex=False)
    report["call"] = (species
                      & (report["top_hit_ani"] >= cfg["skani"]["min_species_ani"])
                      & (report["top_hit_af_query"] >= cfg["skani"]["min_species_af"])
                      & (report["top_hit_af_ref"] >= cfg["skani"]["min_species_af"])
                      ).map({True: "confident_species", False: "not_species"})

    # Each rule is True where the genome FAILS; missing values count as failures
    rules = {
        "missing_seqkit":     report["sum_len"].isna(),
        "length_low":         report["sum_len"] < cfg["seqkit"]["min_length"],
        "length_high":        report["sum_len"] > cfg["seqkit"]["max_length"],
        "n50_low":            report["N50"] < cfg["seqkit"]["min_n50"],
        "missing_checkm2":    report["Completeness"].isna(),
        "completeness_low":   report["Completeness"] < cfg["checkm2"]["min_completeness"],
        "contamination_high": report["Contamination"] > cfg["checkm2"]["max_contamination"],
        "skani_no_hit":       report["top_hit_species"].isna(),
        "not_species":        report["top_hit_species"].notna() & (report["call"] == "not_species"),
    }
    fails = pd.DataFrame(rules)
    report["fail_reasons"] = fails.apply(lambda r: ";".join(fails.columns[r]), axis=1)
    report["qc_pass"] = report["fail_reasons"] == ""

    report.to_csv("qc_report.tsv", sep="\t")
    passed = report[report["qc_pass"]].drop(columns=["qc_pass", "fail_reasons"])
    passed.to_csv("qc_passed.tsv", sep="\t")

    print(f"{len(report)} genomes, {len(passed)} passed, {len(report) - len(passed)} failed")
    print(fails.sum().loc[lambda s: s > 0].to_string())


if __name__ == "__main__":
    main()
