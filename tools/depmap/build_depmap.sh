#!/bin/bash -l
source ~/.bashrc
cd ~/hackatomics/tools/depmap
conda activate jmml_targets

# python build_depmap.py \
#     --data-dir . \
#     --field OncotreePrimaryDisease \
#     --diseases \
#         "Acute Myeloid Leukemia" \
#         "Myeloproliferative Neoplasms" \
#     --subtypes \
#         "Acute Myeloid Leukemia" \
#         "AML, NOS" \
#         "Acute Monoblastic/Monocytic Leukemia" \
#         "AML with Myelodysplasia-Related Changes" \
#         "Myeloid Leukemia Associated with Down Syndrome" \
#     --out depmap_myeloid_features.csv
    
    
    
python build_depmap.py \
    --data-dir . \
    --field OncotreePrimaryDisease \
    --diseases \
        "Fibrosarcoma" \
        "Myxofibrosarcoma" 
    --out depmap_sarcoma_features.csv