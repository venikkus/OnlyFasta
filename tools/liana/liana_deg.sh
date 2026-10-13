#!/bin/bash -l
source ~/.bashrc
cd ~/hackatomics/tools/liana
conda activate jmml_targets

python liana_deg.py --input ../../data/dataset.h5ad --deg ../../data/CD14_Mono_DEGs_surface_GTEx_v11.csv --outdir ../../output/liana_results/CD14_surface/
python liana_deg.py --input ../../data/dataset.h5ad --deg ../../data/HSC_DEGs_surface_GTEx_v11.csv --outdir ../../output/liana_results/HSC_surface/
python liana_deg.py --input ../../data/dataset.h5ad --deg ../../data/GMP_DEGs_surface_GTEx_v11.csv --outdir ../../output/liana_results/GMP_surface/