#!/usr/bin/bash

threCDI=0.5
threEEI=0.5
python ../eeisp.py data.txt Sample --threCDI $threCDI --threEEI $threEEI

python ../gene_name_id.py Sample genes.tsv \
       Sample_CDI_score_data_thre$threCDI.txt \
       Sample_EEI_score_data_thre$threEEI.txt \
       Sample_number_nonzero_exp.txt 0.5
