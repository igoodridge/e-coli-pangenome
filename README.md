# e-coli-pangenome

A Nextflow pipeline that builds an *E. coli* pangenome from genome assemblies, annotates it with COG functional categories, and picks a small set of maximally distant genomes for given COG gene sets.

```
assemblies (*.fasta.gz)  ──┐
extra genomes (e.g. K-12) ─┴─> QC (seqkit, CheckM2, skani) ─> Prokka ─> Panaroo ─> eggNOG-mapper
                                                                                       │
                       genome selection <─ Jaccard distances + PCoA per COG set <─ gene annotation summary
```

1. **QC**: assembly stats (seqkit), completeness/contamination (CheckM2) and species check (skani vs GTDB), filtered with the thresholds in [conf/qc.yml](conf/qc.yml).
2. **Annotation**: Prokka on every genome that passes QC.
3. **Pangenome**: Panaroo (`--clean-mode moderate`).
4. **Functional annotation**: eggNOG-mapper on the Panaroo pangenome reference, with COG IDs described using the NCBI COG 2020 definitions.
5. **Genome selection**: for each COG gene set in [conf/cog_sets.yml](conf/cog_sets.yml), Jaccard distances and a PCoA between genomes, then the genomes that are furthest apart are picked using the rules in [conf/picking.yml](conf/picking.yml).

## Installation

Requirements:

- [Nextflow](https://www.nextflow.io/) ≥ 25.10
- [Conda](https://docs.conda.io/) (or mamba). Nextflow creates each tool's environment from [envs/](envs/) on first use and caches it in `~/.nextflow-conda-cache`.

```bash
git clone <repo-url> e-coli-pangenome
cd e-coli-pangenome
```

### Databases

These are downloaded into `databases/` (set with `--db_cache`) on first run, unless you pass an existing copy:

| Database | Parameter | Notes |
|---|---|---|
| eggNOG (v7.0.0) | `--eggnog_db` | Directory containing `eggnog.db`, `eggnog_proteins.dmnd`, etc. The version is checked before use. |
| CheckM2 | `--checkm2_db` | Path to `uniref100.KO.1.dmnd` |
| skani GTDB r226 | `--skani_db` | Directory of the skani sketch database |

They're large, so pointing to an existing copy is recommended where one exists.

To check the eggNOG environment resolves before a full run:

```bash
conda env create -f envs/eggnog.yml
conda activate eggnog
which emapper.py   # should point inside the eggnog env
```

## Quick example

Run on the first 3 assemblies, with K-12 added as a reference genome:

```bash
nextflow run main.nf \
    --input /path/to/assemblies \
    --extra_genomes assets/extra_genomes.tsv \
    --eggnog_db /path/to/eggnog_db \
    -profile test
```

The `test` profile limits the run to 3 assemblies (`--test_n 3`), lowers CPU use and writes to `results_test/`. Drop `-profile test` for a full run into `results/`.

Add `-resume` to rerun only the steps whose inputs have changed (e.g. editing [conf/picking.yml](conf/picking.yml) only reruns the picking and plots).

## Usage

```bash
nextflow run main.nf --input <dir> [options] [-profile <profile>] [-resume]
```

### Main parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | *required* | Directory of assemblies (`*.fasta.gz`). The file name (minus extension) is used as the genome ID. |
| `--extra_genomes` | none | TSV of additional genomes to include, e.g. K-12 (see [Inputs](#inputs)) |
| `--phylogroups` | none | TSV of phylogroups, used to colour plots and optionally constrain picks |
| `--test_n` | `0` (all) | Use only the first N assemblies |
| `--qc_config` | `conf/qc.yml` | QC thresholds |
| `--cog_sets` | `conf/cog_sets.yml` | COG gene sets to compute distances for |
| `--picking` | `conf/picking.yml` | Genome picking rules |
| `--prokka_dir` | none | Existing Prokka output (`<dir>/<id>/<id>.gff`); genomes found here aren't re-annotated |
| `--qc_seqkit`, `--qc_checkm2`, `--qc_skani` | none | Existing QC reports to reuse instead of rerunning those tools (extra genomes are still run) |
| `--eggnog_db`, `--checkm2_db`, `--skani_db` | none | Existing databases (see [Databases](#databases)) |
| `--db_cache` | `databases/` | Where downloaded databases are stored |
| `-output-dir` | `results` | Where results are published |

The full list, with comments, is in the `params` block of [main.nf](main.nf).

### Profiles

| Profile | Purpose |
|---|---|
| `test` | First 3 assemblies, fewer CPUs, output to `results_test/` |
| `vm` | Input, database and output paths on the analysis VM |
| `vm_reports` | Reuse the QC reports and Prokka output already generated on the VM |

Profiles can be combined, e.g. `-profile vm,vm_reports`.

Resources are set in [nextflow.config](nextflow.config): the local executor is capped at 64 CPUs, with Prokka using 4 CPUs per genome and Panaroo, eggNOG-mapper, CheckM2 and skani each using 64.

## Inputs

**Assemblies** (`--input`): a directory of gzipped FASTA files, one per genome, e.g. `0007d12b-....fasta.gz`.

**Extra genomes** (`--extra_genomes`, optional): a TSV with a header:

```tsv
genome_id	fasta	phylogroup
K-12	/path/to/GCF_000005845.2_ASM584v2_genomic.fna	A
```

`fasta` can be gzipped or not. Extra genomes go through QC but are kept even if they fail (set `extra_genomes.always_keep: false` in [conf/qc.yml](conf/qc.yml) to filter them), and are never limited by `--test_n`.

**Phylogroups** (`--phylogroups`, optional): a TSV mapping genome IDs to phylogroups:

```tsv
guuid	phylogroup
007d12b9-f903-446a-9d66-184fb751cae9	B2
```

### Config files

- [conf/qc.yml](conf/qc.yml): genome length and N50 limits, CheckM2 completeness/contamination, and skani species/ANI/alignment fraction thresholds.
- [conf/cog_sets.yml](conf/cog_sets.yml): named gene sets, each a list of COG category letters. A gene is in a set if any of its categories is listed; `all` means every gene.
- [conf/picking.yml](conf/picking.yml): which COG sets to pick in, a fixed genome to always include (default `K-12`), a minimum gene count, genomes to exclude, and optional one-per-phylogroup constraints.

## Outputs

Everything is published under `results/` (or the `-output-dir` / profile output directory):

```
results/
├── qc/
│   ├── qc_report.tsv          # every genome: metrics, qc_pass, fail_reasons
│   ├── qc_passed.tsv          # genomes that passed QC
│   └── reports/               # seqkit, CheckM2 and skani reports generated in this run
├── prokka/<genome_id>/        # Prokka annotation per genome (.gff, .faa, .gbk, ...)
├── panaroo/                   # Panaroo output, incl. gene_presence_absence.{csv,Rtab},
│                              #   pan_genome_reference.fa, summary_statistics.txt
├── eggnog/                    # eggNOG-mapper output (eggnog_pangenome.emapper.*)
├── annotations/
│   └── gene_annotations.tsv   # one row per gene cluster: frequency, pangenome category
│                              #   (core/soft_core/shell/cloud), eggNOG + COG annotation
├── genome_selection/
│   ├── distances/             # per COG set: <set>.npz (Jaccard distances), <set>.pcoa.tsv,
│   │                          #   plus variance.tsv and sets.tsv
│   ├── picks.tsv              # picked genomes per COG set
│   ├── picks_summary.tsv      # genome length and gene counts of each pick
│   ├── plots/                 # <set>.png PCoA plots and plot_info.tsv
│   └── plots_supp/            # <set>_variance.png, variance explained per component
└── pipeline_info/
    └── software_versions.tsv  # tool versions used
```

## Downstream analysis

Scripts and notebooks outside the pipeline, run on its outputs:

- [scripts/plot_cog_distribution.py](scripts/plot_cog_distribution.py): bar chart of gene clusters per COG category, split into core and accessory.

  ```bash
  python scripts/plot_cog_distribution.py \
      --annotations results/annotations/gene_annotations.tsv \
      --out cog_distribution.png
  ```

- [notebooks/](notebooks/): QC exploration and eggNOG annotation analysis.

## Repository layout

```
e-coli-pangenome/
├── main.nf                          # pipeline entry point and parameters
├── nextflow.config                  # resources, conda envs, profiles
├── modules/                         # one file per pipeline stage
│   ├── 00_qc.nf                     #   seqkit, CheckM2, skani, QC filter
│   ├── 01_annotate.nf               #   Prokka
│   ├── 02_pangenome.nf              #   Panaroo
│   ├── 03_eggnog.nf                 #   eggNOG DB download + eggNOG-mapper
│   ├── 04_summarise.nf              #   merge pangenome with annotations
│   └── 05_genome_selection.nf       #   distances, picking, plots, summary
├── subworkflows/
│   ├── qc.nf                        # runs/reuses QC reports and filters genomes
│   └── genome_selection.nf          # distances → picks → plots → summary
├── bin/                             # Python scripts called by the pipeline
│   ├── qc_filter.py
│   ├── check_eggnog_db.py
│   ├── summarise_annotations.py
│   ├── genome_distances.py
│   ├── pick_genomes.py
│   ├── plot_pcoa.py
│   └── summarise_picks.py
├── conf/
│   ├── qc.yml                       # QC thresholds
│   ├── cog_sets.yml                 # COG gene sets for distances
│   └── picking.yml                  # genome picking rules
├── envs/                            # one conda environment per tool
│   ├── seqkit.yml, checkm2.yml, skani.yml, qc.yml
│   ├── prokka.yml, panaroo.yml, eggnog.yml
│   └── analysis.yml
└── assets/
    ├── cog-20.def.tab               # NCBI COG 2020 definitions
    └── extra_genomes.tsv            # extra genomes sheet (K-12)
```
