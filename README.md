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
How to run EEISP  
`eeisp.sh` performs the calculation of the CDI (Co-Dependency Index)[1] and EEI scores for gene pairs in two steps.  
   1.  `eeisp.py` calculates the CDI and EEI scores for each gene pair. It outputs lists of gene pairs of CDI and EEI, and the tables of degree distribution.     
       `<filename>_CDI_score_data_thre10.0.txt` shows the lists of gene pairs with CDI scores. 
   2.  `gene_name_id.py`converts to the numbers of CDI and EEI gene pairs to Gene Names (Symbols), if the a list of Ensemble Gene IDs 
        and Gene Names are provided. When only Gene Names (Symbols) is provided, `gene_name.py` can be performed.  
```
> sh eeisp.sh  　

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
  * 10.0             # A threhsold for CDI (or EEI)

 ```
### Reference 
[1] S. Mohammadi, J. D. Velderrain, M. Kellis and A. Grama, "DECODE-ing sparsity patterns in single-cell RNA-seq", bioRxiv, 2018. 
