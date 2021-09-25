# EEISP

EEISP identifies gene pairs that are codependent and mutually exclusive from single-cell RNA-seq data. 
       
## Installation
EEISP is written in Python3 and does not require an installation.  

### Dependencies
EEISP requires the following libraries.
* numpy
* pandas
* scipy
* math
* multiprocessing
* time
* cupy (when using GPU computation `--gpu`)

## Usage
EEISP takes a read count matrix as an input, in which rows and columns represent genes and cells, respectively.  

   1.  `eeisp.py` calculates the CDI and EEI scores for all gene pairs. The output contains lists of gene pairs that have CDI or EEI values above the specified threshold and the tables of degree distribution.
       ```
         usage: eeisp [-h] [--threCDI THRECDI] [--threEEI THREEEI] [--tsv] [--gpu] [-p THREADS] [-v] matrix output

         positional arguments:
           matrix                Input matrix
           output                Output prefix

         optional arguments:
           -h, --help            show this help message and exit
           --threCDI THRECDI     Threshold for CDI (default: 20.0)
           --threEEI THREEEI     Threshold for EEI (default: 10.0)
           --tsv                 Specify when the input file is tab-delimited (.tsv)
           --gpu                 GPU mode
           -p THREADS, --threads THREADS  number of threads (default: 2)
           -v, --version         show program's version number and exit
       ```  
   2.  `add_genename_from_geneid.py` add Gene Names (Symbols) to the output files of `eeisp.py`.
        ```
         usage: add_genename_from_geneid.py [-h] [--i_id I_ID] [--i_name I_NAME] input output genelist

         positional arguments:
           input            Input matrix
           output           Output prefix
           genelist         Gene list

         optional arguments:
           -h, --help       show this help message and exit
           --i_id I_ID      column number of gene id (default: 0)
           --i_name I_NAME  column number of gene name (default: 1)
       ```
### Example
The sample data is included in `sample` directory. 
   * `data.txt` is the input matrix of scRNA-seq data.
   * `genelidlist.txt` is the gene list for `add_genename_from_geneid.py`.


```
 eeisp.py data.txt Sample --threCDI 0.5 --threEEI 0.5 -p 8
 add_genename_from_geneid.py Sample_CDI_score_data_thre0.5.txt Sample_CDI_score_data_thre0.5.addgenename.txt geneidlist.txt
 add_genename_from_geneid.py Sample_EEI_score_data_thre0.5.txt Sample_EEI_score_data_thre0.5.addgenename.txt geneidlist.txt
```

* Output files  
```
   Sample_CDI_score_data_thre0.5.txt            # A list of gene pairs with CDI score.  
   Sample_CDI_degree_distribution_thre0.5.csv   # A table of the number of CDI degree and genes.  
   Sample_CDI_score_data_thre0.5.addgenename.txt  # with Gene names
   Sample_EEI_score_data_thre0.5.txt            # A list of gene pairs with EEI scores.  
   Sample_EEI_degree_distribution_thre0.5.csv   # A table of the number of EEI degree and genes.
   Sample_EEI_score_data_thre0.5.addgenename.txt  # with gene names
```

## Reference
Nakajima N., Hayashi T., Fujiki K., Shirahige K., Akiyama T., Akutsu T. and Nakato R., [Codependency and mutual exclusivity for gene community detection from sparse single-cell transcriptome data](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkab601/6324613), *Nucleic Acids Research*, 2021.
