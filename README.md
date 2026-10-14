# JMML Target Evidence Pipeline

This repository contains a configurable Snakemake workflow that starts from an annotated, integrated AnnData object and builds a gene-level evidence matrix for target prioritization. It runs pseudobulk differential expression from raw counts, then joins independent functional, regulatory, cell-surface, target-annotation and compound evidence, and ends with one table of candidate targets. It was developed for JMML, but runs on any single-cell AnnData with a case/control contrast: all dataset-specific settings live in one config file. The output is for research hypothesis generation and experimental prioritization; it does not establish drug efficacy or clinical suitability.

## What is implemented

The workflow in `workflow/Snakefile` orchestrates these stages:

1. **Pseudobulk DE (`tools/differential_expression/pseudobulk_de.py`)** — aggregates raw counts by biological sample and configured cell population, then fits PyDESeq2 with `~ condition` (or `~ subject + condition` when `pseudobulk_de.block_col` defines paired samples). In the JMML example, `PTPN11`, `NF1`, `NRAS`, `KRAS`, and `CBL` are grouped as `JMML`; `Control` is the reference group; `RALD` is excluded from this contrast. The default output filter is `padj < 0.05` and `|log2FC| >= 1`.
2. **DepMap** — intersects DE genes with the configured CRISPR dependency feature table. The JMML example reuses the existing myeloid snapshot; its release metadata is currently unverified.
3. **Surface evidence** — downloads/uses Human Protein Atlas, CSPA, TCSA, and GTEx tables to annotate candidate genes. This is an expression/druggability clue, not proof of tumour specificity or normal-tissue safety.
4. **CollecTRI** — calculates sample-level TF activity from expression profiles and creates candidate TF activity summaries. Condition labels are mapped to case/control through `pseudobulk_de.group_map`.
5. **LIANA** — compares ligand–receptor scores for the configured receiver population between case and control cells, then intersects receptor genes with DE genes. The current implementation pools cells by condition; its score difference is exploratory and is not a patient-level differential test.
6. **Open Targets and ChEMBL** — query target annotation (tractability, genetic constraint, gene name, Reactome pathways) and target–compound evidence for the DE gene list through their APIs.
7. **Integration/report** — joins the evidence by gene without replacing missing values with zero and writes per-population evidence tables plus a combined report.
8. **Summary** — collects candidate genes and the key evidence of every module in one table (see [Outputs](#outputs)).

`CellOracle` is optional and disabled in the included configs. It is not part of the integrated target-evidence table, but its rank and delta score appear in the summary when it runs. The workflow does **not** calculate one summed therapeutic score: the output is an evidence matrix ordered by DE `padj`, so therapeutic ranking still requires transparent interpretation of the separate evidence fields.

## Installation

1. Clone this repository:
```bash
git clone git@github.com:venikkus/OnlyFasta.git
cd OnlyFasta
```

2. Create the runner environment with Snakemake (requires conda or mamba; see the [Miniconda installation instructions](https://www.anaconda.com/docs/getting-started/miniconda/install)). Each workflow step creates its own environment from `workflow/envs/` on first run.
```bash
mamba env create -f workflow/envs/runner.yaml -n onlyfasta
conda activate onlyfasta
```

## Run

Edit a config (see below), preview the jobs with a dry run, then run:

```bash
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --sdm conda --cores 8 -n
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --sdm conda --cores 8
```

On a shared server, limit memory and concurrent API jobs, and run long jobs in `tmux`:

```bash
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --sdm conda \
  --cores 8 --resources mem_mb=40000 api_slots=1 --keep-going --rerun-incomplete
```

On `deeploid`, keep Conda environments, package cache, and temporary files on `/media/hdd` because the system partition is nearly full. For example, add `--conda-prefix /media/hdd/metagenomics/hackatomics/snakemake_envs/OnlyFasta_try` and set `CONDA_PKGS_DIRS` and `TMPDIR` to directories on that mount.

Relative paths resolve from the project root, so Snakemake can be launched from any directory. To run a single module, use the helper scripts in `tools/*/*.sh`; they read `config/config.yaml` by default, or the file given in `CONFIG`:

```bash
CONFIG=config/examples/sarcoma_de.yaml tools/depmap/run_depmap.sh
```

Individual Python tools can also be run directly; pass input and output paths explicitly (see `--help` of each tool).

## Configure for your AnnData

`config/config.yaml` is a commented template; filled-in examples are in `config/examples/` (`jmml_from_adata.yaml`: JMML vs Control per cell type; `sarcoma_de.yaml`: tumour vs microenvironment within samples). For a new dataset, set at least:

| Field | What to set |
|---|---|
| `run_name` | A new name for each analysis; outputs go to `results/<run_name>/` |
| `adata.path`, `adata.counts_layer` | Your `.h5ad` and its raw integer counts (a layer name, or `X`) |
| `adata.sample_column`, `condition_column`, `celltype_column` | Column names in `adata.obs` |
| `adata.case`, `adata.control` | Contrast labels after `pseudobulk_de.group_map` is applied |
| `populations`, `datasets.<pop>.celltype` | Populations to analyse and their cell type labels |
| `pseudobulk_de.group_map`, `cell_types`, `case_label`, `control_label` | Condition mapping and cell types for DE; `block_col` for paired designs |
| `depmap.features` or `depmap.mode: build` + `diseases` | DepMap models relevant to your disease |
| `surface.exclude_gtex_tissues` | GTEx tissues ignored in the normal-tissue filter (blood-forming tissues for haematological tumours, `[]` for solid tumours) |

Each module has a `mode`: `compute` runs it, `reuse` takes a precomputed table from `datasets.<pop>`, and `"off"` (where available) skips it; its columns in the summary stay empty. Keep `"off"` quoted in YAML so it is parsed as a string. Reused resources remain external inputs; no large AnnData or result tables are included in this repository.

### Before you run: checklist

1. **One config = one analysis.** Use a new `run_name` for every analysis; relative paths resolve from the repository root, absolute paths are accepted. Keep machine-specific configs (with your own `adata.path`) outside the repository, e.g. `cp config/examples/sarcoma_de.yaml ~/my_run.yaml`.
2. **Raw integer counts.** `adata.counts_layer` must point to unnormalized integer counts (a key of `adata.layers`, or `X`). Normalized, log or SCT values make the DE step fail. Objects converted from Seurat often keep normalized data in `X`; check before running.
3. **Choose the design; it decides which modules can run.**

   | Design | Example | `pseudobulk_de.group_col` | `block_col` | Modules |
   |---|---|---|---|---|
   | Between patient groups: the same cell type in case vs control patients | JMML: HSC of JMML patients vs HSC of controls | column with the patient group | `null` | all |
   | Within samples: different cell types of the same patients | Sarcoma: tumour cells vs microenvironment | the cell type column | patient / sample ID | switch `collectri`, `liana` and `celloracle` to `"off"`: they compare one cell type between conditions |

4. **`group_map` translates labels to case/control.** Keys are values found in `group_col`, values are the contrast labels (e.g. `PTPN11: JMML`, `Macrophages: Microenvironment`). Values missing from the map are excluded. The resulting labels must equal both `adata.case`/`adata.control` and `pseudobulk_de.case_label`/`control_label`.
5. **Three lists must match:** `populations`, the keys of `datasets` and the keys of `pseudobulk_de.cell_types`. `datasets.<pop>.celltype` is the value of `adata.celltype_column` for that population.
6. **Keep duplicated settings consistent:** `pseudobulk_de.sample_col`, `celltype_col`, `group_col` and `counts_layer` usually equal `adata.sample_column`, `celltype_column`, `condition_column` and `counts_layer`.
7. **Biological replicates are samples, not cells.** Each sample ID must map consistently to one condition (unpaired designs) or to one patient (`block_col`, paired designs); at least `pseudobulk_de.min_replicates` (default 2) samples per condition.
8. **Disease-specific settings:** DepMap models of your disease (`depmap.features`, or `depmap.mode: build` with `diseases`); GTEx tissues ignored in the normal-tissue filter (`surface.exclude_gtex_tissues`: blood-forming tissues for haematological tumours, `[]` for solid tumours); `surface.paper_check: false` unless the data are JMML.
9. **Data requirements:** gene symbols (HGNC) as unique `var_names` (surface, Open Targets and CollecTRI match on symbols); internet access for conda environments, surface atlases, Open Targets and ChEMBL.
10. **Resources:** `threads` and `mem_mb` per module are hints for the scheduler. On a laptop run without `--resources`; on a server, `--resources mem_mb=<available>` keeps parallel jobs within memory.
11. **Always dry-run first** (`-n`), then run; read `results/<run_name>/summary/` first.

To check a dataset before writing the config, list its metadata and test the counts without loading the whole matrix (run in the `onlyfasta-deseq2` environment: `mamba env create -f workflow/envs/deseq2.yaml -n onlyfasta-deseq2`):

```bash
python - <<'EOF'
import sys, numpy as np
sys.path.insert(0, "tools/common")
from lazy_h5ad import open_counts
handle, adata = open_counts("path/to/your_dataset.h5ad", "X")   # or the name of the counts layer
print(adata.n_obs, "cells,", adata.n_vars, "genes; obs columns:", list(adata.obs.columns))
rows = (adata.X if adata.X is not None else next(iter(adata.layers.values())))[np.arange(min(500, adata.n_obs)), :]
values = rows.data if hasattr(rows, "data") else np.asarray(rows).ravel()
print("raw integer counts:", bool(np.all(values >= 0) and np.all(values == np.round(values))))
handle.close()
EOF
```

### Input contract of the JMML example

`config/examples/jmml_from_adata.yaml` expects an annotated `.h5ad` with:

- `adata.layers["counts"]`: non-negative integer raw counts;
- `adata.obs["batch"]`: one independent patient/sample ID per biological replicate;
- `adata.obs["sample_groups"]`: mutation/control labels;
- `adata.obs["celltype_l2_merged"]`: cell annotations.

The included populations are `HSC`, `GMP`, and `Mono_CD14`, mapped to the annotation values `HSC`, `GMP`, and `CD14 Mono`. The workflow uses `sample_groups` because the alternative `sample_groups_wide` field is incomplete for some control samples in this atlas. All controls mapped to `Control` are included in the primary contrast; if controls 10–12 are considered atypical, run a separate sensitivity analysis and compare it with the primary result. Each `batch` must map consistently to one condition.

## Outputs

Everything is written to `results/<run_name>/` (under `paths.results_root`). **Start with `summary/`**: it collects the candidates and the key evidence of every module in one table.

- `summary/final_candidates.csv` and `.xlsx` (sheet `all` + one per population): one row per candidate gene. Candidates are the union of the DEG branch (upregulated, membrane, low normal-tissue expression) and the regulatory branch (DEG transcription factors with higher CollecTRI activity, CellOracle hubs). Columns: gene name and Reactome pathways, DE log2FC and adjusted p-value, surface atlases, tumor specificity and GTEx TPM, top-5 LIANA interactions, DepMap selectivity, Open Targets, CellOracle rank and delta score, TF activity and p-value.
- `summary/README.md`: description of every column.
- `summary/opentargets_report_examples.md`: Open Targets report for the top candidates (`summary.opentargets_report_genes`).

The summary can also be rebuilt from an existing run: `python tools/summary/build_summary.py --results-dir results/<run_name> --populations HSC GMP --case JMML --control Control --outdir results/<run_name>/summary`.

Per-module outputs stay available for details; each step directory also contains a `manifest.json` with parameters, input checksums, package versions, QC counts and warnings. Stage logs and runtime/resource records go to `logs/<run_name>/` and `benchmarks/<run_name>/` (under `paths.logs_root` and `paths.benchmarks_root`); the Snakemake log of each run is in `.snakemake/log/`.

```
results/<run_name>/
├── summary/                              # final candidate table (see above)
├── pseudobulk_de/<pop>_DEGs.csv          # DE results (pseudobulk_de.mode: compute)
├── de/<pop>/genes.csv                    # normalized DEG table used by all modules
├── depmap/<pop>/features.csv             # CRISPR dependency per gene (+ features_candidates.csv)
├── surface/<pop>/                        # DEGs_surface.csv, DEGs_surface_GTEx.csv, candidates.csv, DEGs_surface_GTEx.xlsx
├── collectri/<pop>/candidate_tf_activity.csv
├── liana/<pop>/JMML_vs_Control_DEG_receptors.csv
├── opentargets/<pop>/opentargets_features.csv
├── chembl/<pop>/chembl_gene_features.csv
├── celloracle/<pop>_rank_candidates.csv  # + celloracle_graph_metrics.csv (not in the integrated table)
├── integrated/<pop>/target_evidence.csv  # one row per DEG with all evidence columns
└── report/                               # all_target_evidence.csv, run_summary.csv, REPORT.md, config.resolved.json
```

Missing evidence remains missing; it is not interpreted as a negative result. A gene appearing in several evidence sources is a stronger hypothesis, but expression, cell-line dependency, ligand–receptor association, target engagement, and clinical activity are different kinds of evidence and should be reviewed separately. Validate promising candidates experimentally.

Example outputs of the JMML analysis are in `examples/`.
