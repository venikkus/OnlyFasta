#!/bin/bash -l
source ~/.bashrc

cd ~/hackatomics/tools/opentarget
conda activate jmml_targets

python parse_opentarget.py --input ../../data/CD14_Mono_candidates_GTEx_v11.csv --report --outdir ../../output/opentargets_results/opentarget_cd14_reports_final --report-max-genes 70
python parse_opentarget.py --input ../../data/HSC_candidates_GTEx_v11.csv --report --outdir ../../output/opentargets_results/opentarget_hsc_reports_final --report-max-genes 70
python parse_opentarget.py --input ../../data/GMP_candidates_GTEx_v11.csv --report --outdir ../../output/opentargets_results/opentarget_gmp_reports_final --report-max-genes 70