# JMML Target Evidence Pipeline (JTEP)

Snakemake workflow for prioritizing candidate therapeutic targets from cell type specific differential expression and supporting molecular evidence. It was developed for JMML, but runs on any single-cell AnnData with a case/control contrast: all dataset-specific settings live in one config file.

JTEP integrates cell type specific differential expression with evidence on tumor selectivity, CRISPR dependency, transcriptional regulation, cell-cell communication, disease associations, and target-drug interactions. It produces traceable evidence tables, ranked candidate targets, and a summary report. Results are hypothesis-generating and require experimental validation; they do not establish drug efficacy or clinical suitability.

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

Edit `config/config.yaml` (see below), preview the jobs with a dry run, then run:

```bash
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --sdm conda --cores 8 -n
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --sdm conda --cores 8
```

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

Each module has a `mode`: `compute` runs it, `reuse` takes a precomputed table from `datasets.<pop>`, and `"off"` (where available) skips it; its columns in the summary stay empty.

### Before you run: checklist

1. **One config = one analysis.** Use a new `run_name` for every analysis; relative paths resolve from the repository root. Keep machine-specific configs (with your own `adata.path`) outside the repository, e.g. `cp config/examples/sarcoma_de.yaml ~/my_run.yaml`.
2. **Raw integer counts.** `adata.counts_layer` must point to unnormalized integer counts (a key of `adata.layers`, or `X`). Normalized, log or SCT values make the DE step fail. Objects converted from Seurat often keep normalized data in `X`; check before running.
3. **Choose the design; it decides which modules can run.**

   | Design | Example | `pseudobulk_de.group_col` | `block_col` | Modules |
   |---|---|---|---|---|
   | Between patient groups: the same cell type in case vs control patients | JMML: HSC of JMML patients vs HSC of controls | column with the patient group | `null` | all |
   | Within samples: different cell types of the same patients | Sarcoma: tumour cells vs microenvironment | the cell type column | patient / sample ID | switch `collectri`, `liana` and `celloracle` to `"off"`: they compare one cell type between conditions |

4. **`group_map` translates labels to case/control.** Keys are values found in `group_col`, values are the contrast labels (e.g. `PTPN11: JMML`, `Macrophages: Microenvironment`). Values missing from the map are excluded. The resulting labels must equal both `adata.case`/`adata.control` and `pseudobulk_de.case_label`/`control_label`.
5. **Three lists must match:** `populations`, the keys of `datasets` and the keys of `pseudobulk_de.cell_types`. `datasets.<pop>.celltype` is the value of `adata.celltype_column` for that population.
6. **Keep duplicated settings consistent:** `pseudobulk_de.sample_col`, `celltype_col`, `group_col` and `counts_layer` usually equal `adata.sample_column`, `celltype_column`, `condition_column` and `counts_layer`.
7. **Disease-specific settings:** DepMap models of your disease (`depmap.features`, or `depmap.mode: build` with `diseases`); GTEx tissues ignored in the normal-tissue filter (`surface.exclude_gtex_tissues`: blood-forming tissues for haematological tumours, `[]` for solid tumours); `surface.paper_check: false` unless the data are JMML.
8. **Data requirements:** gene symbols (HGNC) as unique `var_names` (surface, Open Targets and CollecTRI match on symbols); at least `pseudobulk_de.min_replicates` (default 2) samples per condition; internet access for conda environments, surface atlases, Open Targets and ChEMBL.
9. **Resources:** `threads` and `mem_mb` per module are hints for the scheduler. On a laptop run without `--resources`; on a server, `--resources mem_mb=<available>` keeps parallel jobs within memory.
10. **Always dry-run first** (`-n`), then run; read `results/<run_name>/summary/` first.

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

## Outputs

Everything is written to `results/<run_name>/`. **Start with `summary/`**: it collects the candidates and the key evidence of every module in one table.

- `summary/final_candidates.csv` and `.xlsx` (sheet `all` + one per population): one row per candidate gene. Candidates are the union of the DEG branch (upregulated, membrane, low normal-tissue expression) and the regulatory branch (DEG transcription factors with higher CollecTRI activity, CellOracle hubs). Columns: gene name and Reactome pathways, DE log2FC and adjusted p-value, surface atlases, tumor specificity and GTEx TPM, top-5 LIANA interactions, DepMap selectivity, Open Targets, CellOracle rank and delta score, TF activity and p-value.
- `summary/README.md`: description of every column.
- `summary/opentargets_report_examples.md`: Open Targets report for the top candidates (`summary.opentargets_report_genes`).

The summary can also be rebuilt from an existing run: `python tools/summary/build_summary.py --results-dir results/<run_name> --populations HSC GMP --case JMML --control Control --outdir results/<run_name>/summary`.

Per-module outputs stay available for details; each step directory also contains a `manifest.json` with parameters, input checksums, package versions, QC counts and warnings. Logs and benchmarks go to `logs/<run_name>/` and `benchmarks/<run_name>/`.

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

Example outputs of the JMML analysis are in `examples/`.
