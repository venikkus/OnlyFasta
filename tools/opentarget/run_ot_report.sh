#!/bin/bash -l
source ~/.bashrc

cd ~/hackatomics/tools/opentarget
conda activate jmml_targets

python parse_opentarget.py --input ../../data/CD14_Mono_candidates_GTEx_v11.csv --report --outdir ../../output/opentargets_results/opentarget_cd14_reports --report-max-genes 70