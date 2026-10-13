#!/bin/bash -l
source ~/.bashrc
cd ~/hackatomics/tools/liana
conda activate jmml_targets

python top5_receptors.py \
    --input ../../output/liana_results/CD14_surface/JMML_vs_Control_DEG_receptors.csv \
    --out ../../output/liana_results/CD14_surface/top5_receptors.csv
python top5_receptors.py \
    --input ../../output/liana_results/GMP_surface/JMML_vs_Control_DEG_receptors.csv \
    --out ../../output/liana_results/GMP_surface/top5_receptors.csv
python top5_receptors.py \
    --input ../../output/liana_results/HSC_surface/JMML_vs_Control_DEG_receptors.csv \
    --out ../../output/liana_results/HSC_surface/top5_receptors.csv