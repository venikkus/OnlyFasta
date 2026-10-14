# JMML Target Evidence Pipeline

This repository contains a configurable Snakemake workflow that starts from an annotated, integrated AnnData object and builds a gene-level evidence matrix for target prioritization in JMML. It runs pseudobulk differential expression from raw counts, then joins independent functional, regulatory, cell-surface, disease-association, and compound evidence. The output is for research hypothesis generation and experimental prioritization.

## What is implemented

The workflow in `workflow/Snakefile` orchestrates these stages:

1. **Pseudobulk DE (`tools/differential_expression/pseudobulk_de.py`)** — aggregates raw counts by biological sample and configured cell population, then fits PyDESeq2 with `~ condition`. In the included JMML config, `PTPN11`, `NF1`, `NRAS`, `KRAS`, and `CBL` are grouped as `JMML`; `Control` is the reference group; `RALD` is excluded from this contrast. The default output filter is `padj < 0.05` and `|log2FC| >= 1`.
2. **DepMap** — intersects DE genes with the configured CRISPR dependency feature table. The JMML config reuses the existing myeloid snapshot; its release metadata is currently unverified.
3. **Surface evidence** — downloads/uses Human Protein Atlas, CSPA, TCSA, and GTEx tables to annotate candidate genes. This is an expression/druggability clue, not proof of leukemia specificity or normal-HSPC safety.
4. **CollecTRI** — calculates sample-level TF activity from expression profiles and creates candidate TF activity summaries. The JMML mutation labels are mapped to the common `JMML` group through configuration.
5. **LIANA** — compares ligand–receptor scores for the configured receiver population between JMML and control cells, then intersects receptor genes with DE genes. The current implementation pools cells by condition; its score difference is exploratory and is not a patient-level differential test.
6. **Open Targets and ChEMBL** — queries disease-association and target–compound evidence for the DE gene list through their APIs.
7. **Integration/report** — joins the evidence by gene without replacing missing values with zero and writes per-population evidence tables plus a combined report.

`CellOracle` is optional and disabled in the included config. It is not currently part of the integrated target-evidence table. The workflow does **not** calculate one summed therapeutic score: the output is an evidence matrix ordered by DE `padj`, so therapeutic ranking still requires transparent interpretation of the separate evidence fields.

## Input contract

The input is an annotated `.h5ad`. The JMML configuration expects:

- `adata.layers["counts"]`: non-negative integer raw counts;
- `adata.obs["batch"]`: one independent patient/sample ID per biological replicate;
- `adata.obs["sample_groups"]`: mutation/control labels;
- `adata.obs["celltype_l2_merged"]`: cell annotations.

The included populations are `HSC`, `GMP`, and `Mono_CD14`, mapped to the annotation values `HSC`, `GMP`, and `CD14 Mono`. The workflow uses `sample_groups` because the alternative `sample_groups_wide` field is incomplete for some control samples in this atlas. All controls mapped to `Control` are included in the primary contrast; if controls 10–12 are considered atypical, run a separate sensitivity analysis and compare it with the primary result. Each `batch` must map consistently to one condition. Samples are the biological replicates; cells are not treated as independent replicates.

## Configuration

The dataset-specific example is `config/jmml_from_adata.yaml`. Change `adata.path`, metadata column names, `populations`, `datasets`, and `pseudobulk_de.group_map` / `cell_types` for another dataset. Relative paths resolve from the project root; absolute paths are accepted. `pseudobulk_de.mode: compute` means DE is generated from AnnData during this run, with no precomputed DE files required.

Other stages can be set to `compute`, `reuse`, or `off` where supported. Reused resources must be supplied in the config and remain external inputs; no large AnnData or result tables are included in this repository. Keep `celloracle.mode` quoted as `"off"` in YAML so it is parsed as a string.

## Run

From the repository root, first preview the DAG:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config/jmml_from_adata.yaml \
  --cores 8 --dry-run
```

Run with the workflow-specific Conda environments:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config/jmml_from_adata.yaml \
  --cores 8 --resources mem_mb=40000 api_slots=1 \
  --use-conda
```

On `deeploid`, keep Conda environments, package cache, and temporary files on `/media/hdd` because the system partition is nearly full. For example, add `--conda-prefix /media/hdd/metagenomics/hackatomics/snakemake_envs/OnlyFasta_try` and set `CONDA_PKGS_DIRS` and `TMPDIR` to directories on that mount. Run long jobs in `tmux` and inspect `logs/snakemake/<run_name>/run.log`.

## Outputs

For the included run name, outputs are under `results/jmml_from_adata_20261014/`:

- `pseudobulk_de/{HSC,GMP,Mono_CD14}_DEGs.csv` — DE tables computed from the AnnData counts layer;
- `integrated/<population>/target_evidence.csv` — evidence matrix for each population;
- `report/all_target_evidence.csv` — combined gene-by-population evidence matrix;
- `report/run_summary.csv`, `report/REPORT.md`, and `report/config.resolved.json` — run summary, warnings, and resolved settings;
- `logs/snakemake/<run_name>/` and `benchmarks/snakemake/<run_name>/` — stage logs and runtime/resource records.

Missing evidence remains missing; it is not interpreted as a negative result. A gene appearing in several evidence sources is a stronger hypothesis, but expression, cell-line dependency, ligand–receptor association, target engagement, and clinical activity are different kinds of evidence and should be reviewed separately. Validate promising candidates experimentally.
