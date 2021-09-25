#!/usr/bin/bash

threCDI=0.5
threEEI=0.5
../eeisp.py data.txt Sample --threCDI $threCDI --threEEI $threEEI

../add_genename_from_geneid.py Sample_CDI_score_data_thre$threCDI.txt Sample_CDI_score_data_thre$threCDI.addgenename.txt geneidlist.txt
../add_genename_from_geneid.py Sample_EEI_score_data_thre$threEEI.txt Sample_EEI_score_data_thre$threEEI.addgenename.txt geneidlist.txt
