# JMML Target Evidence Pipeline (JTEP)

Snakemake workflow for prioritizing candidate therapeutic targets in JMML from cell type specific differential expression and supporting molecular evidence.

JTEP integrates cell type specific differential expression with evidence on tumor selectivity, CRISPR dependency, transcriptional regulation, cell-cell communication, disease associations, and target-drug interactions. It produces traceable evidence tables, ranked candidate targets, and a summary report. Results are hypothesis-generating and require experimental validation; they do not establish drug efficacy or clinical suitability.

## Installation

1. Clone this repository to your local machine:
```bash
git clone git@github.com:venikkus/OnlyFasta.git
cd OnlyFasta
```

2. Use provided `environment.yaml` file to create the conda environment. If you don't have conda, use these [installation instuctions](https://www.anaconda.com/docs/getting-started/miniconda/install#linux):

```bash
conda env create -f environment.yaml    
conda activate target_search
```

## Configure paths and run

Edit `config/config.yaml` for a dataset or machine. Relative paths resolve from the `OnlyFasta` project root, so the same config works regardless of the directory from which Snakemake is launched. Set `paths.project_root` only when the project root is outside its default location; set `paths.results_root`, `paths.logs_root`, `paths.benchmarks_root`, and `paths.cache_root` to relocate generated files. `paths.tools` and `paths.workflow` let you relocate scripts and environment definitions. Dataset inputs, including `adata.path`, are configured in the same file.

Preview the configured workflow before running it:

```bash
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --cores 4 --dry-run
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --cores 4
```

The helper scripts under `tools/*/*.sh` locate the project from their own location and run through this config. Direct calls to individual Python tools must pass their input and output paths explicitly; the Snakemake adapters supply those paths from the config.
