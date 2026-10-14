# JMML surface_gtex example

Output of `tools/surface_gtex/surface_gtex_pipeline.py` for JMML vs Control DEGs in HSC, GMP and CD14 Mono
(inputs: `data/HSC_DEGs_Masha.csv`, `data/GMP_DEGs_Masha.csv`, `data/CD14_Mono_DEGs_Masha.csv`).

Parameters: HPA + CSPA + TCSA surface atlases, GTEx v11 median TPM < 20 in normal tissues,
excluding whole blood, spleen and EBV-transformed lymphocytes.

- `DEGs_surface_GTEx_v11.xlsx`: one `<population>_all` sheet (every DEG with atlas and GTEx annotation)
  and one `<population>_candidates` sheet (upregulated, membrane, GTEx-specific) per population.
- The same tables as CSV: `data/<population>_DEGs_surface_GTEx_v11.csv` and
  `data/<population>_candidates_GTEx_v11.csv`.

This example was produced by the pre-workflow version of the tool, which wrote all populations into one
workbook. The workflow now writes per population: `results/<run_name>/surface/<population>/`.
