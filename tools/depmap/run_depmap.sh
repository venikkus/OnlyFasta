#!/bin/bash -l
source ~/.bashrc
cd ~/hackatomics/tools/depmap
conda activate jmml_targets

python run_depmap.py --deg ../../data/CD14_Mono_candidates_GTEx_v11.csv --depmap depmap_features.csv --out ../../output/depmap_results/depmap_results_cd14_candidates/depmap_results_cd14_atlases.csv

python run_depmap.py --deg ../../data/GMP_candidates_GTEx_v11.csv --depmap depmap_features.csv --out ../../output/depmap_results/depmap_results_gmp_candidates/depmap_results_gmp_atlases.csv

python run_depmap.py --deg ../../data/HSC_candidates_GTEx_v11.csv --depmap depmap_features.csv --out ../../output/depmap_results/depmap_results_hsc_candidates/depmap_results_hsc_atlases.csv