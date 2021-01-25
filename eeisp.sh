#!/usr/bin/bash

python3 input_file.csv filename --threCDI 10.0 --threEEI 10.0
python gene_name_id.py filename genes.tsv <filename>_CDI_score_data_thre10.0.txt <filename>_EEI_score_data_thre10.0.txt <filename>_number_nonzero_exp.txt 10.0 
