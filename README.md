# estimateEEI

## Description
This is a tool for identifying co-dependency and mutual exclusivity from sparse single-cell 
RNA-seq data. 


## Usage
### Run estimateEEI
estimateEEI computes gene pairs that are co-dependent and mutually exclusive from 
scRNA-seq data and takes an input as a read count matrix data with rows representing genes and 
columns representing cells. 

### Example
1. How to run estimateEEI  
```
> python estimate_cdi_eei_ver1.py filtered_gene_bc_matrices_P10.csv Stem_P10 --threCDI 10.0 --threEEI 10.0　　

Option
  * filtered_gene_bc_matrices_P10.csv  　　# An input file (.csv) which forms a comma delimited file
  * Stem_P10          # A file name of run
  * --threCDI 10.0    # Set a threshold for CDI(Co-Dependency Index) [1] which users determine
  * --threEEI 10.0    # Set a threshold for EEI which users determine
 ```

2. How to convert numbers to Gene Names (Symbols) 
   1. When Ensemble Gene IDs and Gene Names (Symbols) are given
```
> python gene_name_id.py Stem_P10 genes.tsv Stem_P10_CDI_score_data_thre10.0.txt Stem_P10_number_nonzero_exp.txt 0 10.0
 
Option
  * Stem_P10         # A file name of run
  * genes.tsv        # A list of number, Ensemble Gene IDs and Gene Names (or Symbols), which forms a tab delimited file 
  * Stem_P10_CDI_score_data_thre10.0.txt      # A list of gene pairs with CDI scores
  * Stem_P10_number_nonzero_exp.txt      　　　# A list of genes that have non-zero expressions in scRNA-seq data
  * 0                # A type of a conversion file with CDI (0) or EEI (1)
  * 10.0             # A threhsold for CDI (or EEI)
```
2. 
   2. When only Gene Names (Symbols) is given 
```
> python gene_name.py Stem_P10 genes.tsv Stem_P10_CDI_score_data_thre10.0.txt Stem_P10_number_nonzero_exp.txt 0 10.0

Option
  * genes.tsv        # A list of number and Gene Names (or Symbols), which forms a tab delimited file 
```




### Reference 
[1] S. Mohammadi, J. D. Velderrain, M. Kellis and A. Grama, "DECODE-ing sparsity patterns in single-cell RNA-seq", bioRxiv, 2018. 
