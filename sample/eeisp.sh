#!/usr/bin/bash

python3 eeisp.py data.txt Sample --threCDI 0.5 --threEEI 0.5
python gene_name_id.py Sample genes.tsv Sample_CDI_score_data_thre0.5.txt Sample_EEI_score_data_thre0.5.txt Sample_number_nonzero_exp.txt 0.5

