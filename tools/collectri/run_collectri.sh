#!/bin/bash -l
source ~/.bashrc
cd ~/hackatomics/tools/collectri
conda activate jmml_targets

#python run_collectri.py --input ../../data/dataset.h5ad --outdir ../../data/collectri_data  --min-cells 20

python run_collectri.py --reuse --reuse-dir ../../data/collectri_data/ --outdir ../../output/collectri_results/collectri_hsc_deg_candidates --deg ../../data/HSC_DEGs_Masha.csv --plot-celltype "hsc" --plot-top 5 --network-targets 15 --deg-celltype "hsc" --gene-col 'Unnamed: 0'
python run_collectri.py --reuse --reuse-dir ../../data/collectri_data/ --outdir ../../output/collectri_results/collectri_cd14_deg_candidates --deg ../../data/CD14_Mono_DEGs_Masha.csv --plot-celltype "cd14 mono" --plot-top 5 --network-targets 15 --deg-celltype "cd14 mono" --gene-col 'Unnamed: 0'
python run_collectri.py --reuse --reuse-dir ../../data/collectri_data/ --outdir ../../output/collectri_results/collectri_gmp_deg_candidates --deg ../../data/GMP_DEGs_Masha.csv --plot-celltype "gmp" --network-targets 15 --deg-celltype "gmp" --gene-col 'Unnamed: 0' --plot-top 5