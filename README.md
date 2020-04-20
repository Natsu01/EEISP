# EEISP

## Description
This is a tool for identifying codependent and mutually exclusive gene sets from sparse single-cell 
RNA-seq data. 


## Usage
### Run EEISP
EEISP identifies gene pairs that are codependent and mutually exclusive from sparse 
scRNA-seq data and takes as input a read count matrix data with rows representing genes and 
columns representing cells. 

### Example
1. How to run EEISP  
```
> python eeisp.py filtered_gene_bc_matrices_P10.csv Stem --threCDI 10.0 --threEEI 10.0　　

Option
  * filtered_gene_bc_matrices_P10.csv  　　# An input file forms a comma delimited file (.csv)
  * Stem          # A file name of run
  * --threCDI 10.0    # Set a threshold for CDI(Co-Dependency Index) [1] which users determine
  * --threEEI 10.0    # Set a threshold for EEI which users determine
 ```

2. How to convert numbers to Gene Names (Symbols) 
   1. When Ensemble Gene IDs and Gene Names (Symbols) are provided
```
> python gene_name_id.py Stem genes.tsv Stem_P10_CDI_score_data_thre10.0.txt Stem_P10_number_nonzero_exp.txt 0 10.0
 
Option
  * Stem         # A file name of run
  * genes.tsv        # A list of numbers, Ensemble Gene IDs and Gene Names (or Symbols), which forms a tab delimited file 
  * Stem_P10_CDI_score_data_thre10.0.txt      # A list of gene pairs with CDI scores
  * Stem_P10_number_nonzero_exp.txt      　　　# A list of genes that show nonzero expression in scRNA-seq data
  * 0                # A type of a conversion file with CDI (0) or EEI (1)
  * 10.0             # A threhsold for CDI (or EEI)
```
2. 
   2. When only Gene Names (Symbols) is provided 
```
> python gene_name.py Stem genes.tsv Stem_P10_CDI_score_data_thre10.0.txt Stem_P10_number_nonzero_exp.txt 0 10.0

Option
  * genes.tsv        # A list of numbers and Gene Names (or Symbols), which forms a tab delimited file 
```




### Reference 
[1] S. Mohammadi, J. D. Velderrain, M. Kellis and A. Grama, "DECODE-ing sparsity patterns in single-cell RNA-seq", bioRxiv, 2018. 
