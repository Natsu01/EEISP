#!/usr/bin/bash

threCDI=0.5
threEEI=0.5
eeisp data.txt Sample --threCDI $threCDI --threEEI $threEEI

eeisp_add_genename_from_geneid Sample_CDI_score_data_thre$threCDI.txt Sample_CDI_score_data_thre$threCDI.addgenename.txt geneidlist.txt
eeisp_add_genename_from_geneid Sample_EEI_score_data_thre$threEEI.txt Sample_EEI_score_data_thre$threEEI.addgenename.txt geneidlist.txt
