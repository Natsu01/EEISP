# EEISP

## Description
This is a tool for identifying codependent and mutually exclusive gene sets from sparse single-cell 
RNA-seq data. 


## Usage
EEISP identifies gene pairs that are codependent and mutually exclusive from sparse 
scRNA-seq data and takes as input a read count matrix data with rows representing genes and 
columns representing cells.  

How to run EEISP  
`eeisp.sh` performs the calculation of the CDI (Co-Dependency Index) and EEI scores for gene pairs in two steps.  
   1.  `eeisp.py` calculates the CDI and EEI scores for each gene pair. It outputs lists of gene pairs of CDI and EEI, and the tables of degree distribution.    
       ```
         * input_file.csv    # An input file forms a comma delimited file (.csv).
         * filename          # A file name of run.
         * --threCDI 10.0    # Set a threshold for CDI which users determine.
         * --threEEI 10.0    # Set a threshold for EEI which users determine.
       ```  
   2.  `gene_name_id.py`converts to the numbers of CDI and EEI gene pairs to Gene Names (Symbols), if the a list of Ensemble Gene IDs 
        and Gene Names is provided. When only Gene Names (Symbols) is provided, `gene_name.py` can be performed.  
        ```
         * genes.tsv        # A list of numbers, Ensemble Gene IDs and Gene Names (or Symbols), which forms a tab delimited file. 
         * <filename>_CDI_score_data_thre10.0.txt      # A list of gene pairs with CDI scores.  
         * <filename>_EEI_score_data_thre10.0.txt      # A list of gene pairs with EEI scores. 
         * <filename>_number_nonzero_exp.txt      　　　# A list of genes that are expressed in more than at least one cell.
         * 10.0             # A threhsold for CDI (or EEI).
       ```
### Installation
EEISP is written in Python3 and can run without instllation.  


### Example
The sample data is included in `sample`. 
   * `data.txt` The input matrix of scRNA-seq data.
```
 ## eeisp.sh
 #!/usr/bin/bash

 python3 eeisp.py data.txt Sample --threCDI 0.5 --threEEI 0.5
 python gene_name_id.py Sample genes.tsv Sample_CDI_score_data_thre0.5.txt Sample_EEI_score_data_thre0.5.txt Sample_number_nonzero_exp.txt 0.5
```
```
 ## To perform the calculation of CDI and EEI. 
> sh eeisp.sh 
```
  * Output files  
```
   Sample_CDI_score_data_thre0.5.txt            # A list of gene pairs with CDI score.  
   Sample_CDI_degree_distribution_thre0.5.csv   # A table of the number of CDI degree and genes.  
   Sample_CDI_convert_data_thre0.5.txt          # A converted file of the CDI score data.  
   Sample_EEI_score_data_thre0.5.txt            # A list of gene pairs with EEI scores.  
   Sample_EEI_degree_distribution_thre0.5.csv   # A table of the number of EEI degree and genes.   
   Sample_EEI_convert_data_thre0.5.txt          # A converted file of the EEI score data.  
```
<!--###
Option
* eeisp.py
  * input_file.csv    # An input file forms a comma delimited file (.csv).
  * filename          # A file name of run.
  * --threCDI 10.0    # Set a threshold for CDI which users determine.
  * --threEEI 10.0    # Set a threshold for EEI which users determine.
 
* gene_name_id.py 
  * genes.tsv        # A list of numbers, Ensemble Gene IDs and Gene Names (or Symbols), which forms a tab delimited file. 
  * <filename>_CDI_score_data_thre10.0.txt      # A list of gene pairs with CDI scores.  
  * <filename>_EEI_score_data_thre10.0.txt      # A list of gene pairs with EEI scores. 
  * <filename>_number_nonzero_exp.txt      　　　# A list of genes that are expressed in more than at least one cell.
  * 10.0             # A threhsold for CDI (or EEI).

 ```--->
<!--### Reference)
[]([1] S. Mohammadi, J. D. Velderrain, M. Kellis and A. Grama, "DECODE-ing sparsity patterns in single-cell RNA-seq", bioRxiv, 2018.--->
