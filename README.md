# EEISP

EEISP identifies gene pairs that are codependent and mutually exclusive from single-cell RNA-seq data. 
       
## Installation

    pip install eeisp

###  (Optional) Dependencies for GPU
EEISP requires [cupy](https://cupy.dev/) when using GPU computation `--gpu`. Use pip to install cupy like this (see [the manual](https://docs.cupy.dev/en/stable/install.html) for more detail).

    # For CUDA 9.2
    pip install cupy-cuda92
    # For CUDA 10.1
    pip install cupy-cuda101

If you do not use `--gpu`, you do not need to install cupy.

## Usage
EEISP takes a read count matrix as an input, in which rows and columns represent genes and cells, respectively.  

   0. (Optional) Convert CellRanger output to an input matrix (require R and [Seurat](https://satijalab.org/seurat/) library)
       ```
         datadir="outs/filtered_feature_bc_matrix/"
         matrix="matrix.txt"
         R -e "library(Seurat); so <- Read10X('$datadir'); write.table(so, '$matrix', quote=F, sep=',', col.names=T)"
       ```

   1.  `eeisp` calculates the CDI and EEI scores for all gene pairs. The output contains lists of gene pairs that have CDI or EEI values above the specified threshold and the tables of degree distribution.
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
   2.  `eeisp_add_genename_from_geneid` add Gene Names (Symbols) to the output files of `eeisp`.
        ```
         usage: eeisp_add_genename_from_geneid [-h] [--i_id I_ID] [--i_name I_NAME] input output genelist

         positional arguments:
           input            Input matrix
           output           Output prefix
           genelist         Gene list

         optional arguments:
           -h, --help       show this help message and exit
           --i_id I_ID      column number of gene id (default: 0)
           --i_name I_NAME  column number of gene name (default: 1)
       ```
## Tutorial
The sample data is included in `sample` directory. 
   * `data.txt`: the input matrix of scRNA-seq data.
   * `genelidlist.txt`: the gene list for `eeisp_add_genename_from_geneid`.


    eeisp data.txt Sample --threCDI 0.5 --threEEI 0.5 -p 8
This command outputs gene pair lists that have CDI>0.5 or EEI>0.5. `-p 8` means 8 CPUs are used.

Supply `--gpu` option to GPU computation (require [cupy](https://www.preferred.jp/en/projects/cupy/)):

    eeisp data.txt Sample --threCDI 0.5 --threEEI 0.5 -p 8 --gpu
    
(Note: Since GPU computation covers a part of eeisp, it is better to use multiple CPUs even in `--gpu` mode for the fast computation.)

Output files are:
```
   Sample_CDI_score_data_thre0.5.txt            # A list of gene pairs with CDI score.  
   Sample_CDI_degree_distribution_thre0.5.csv   # A table of the number of CDI degree and genes.  
   Sample_EEI_score_data_thre0.5.txt            # A list of gene pairs with EEI scores.  
   Sample_EEI_degree_distribution_thre0.5.csv   # A table of the number of EEI degree and genes.
```
The output files include gene ids only.

```
   $ head Sample_CDI_score_data_thre0.5.txt
   2       7       ESG000003       ESG000008       0.96384320244841
   0       1       ESG000001       ESG000002       0.6852891560232545
   0       6       ESG000001       ESG000007       0.6852891560232545
   7       8       ESG000008       ESG000009       0.6852891560232545
   3       9       ESG000004       ESG000010       0.6469554204484568
   4       6       ESG100005       ESG000007       0.5258703930217091
```

Use `eeisp_add_genename_from_geneid` to add gene names using `geneidlist.txt`, which contains the pairs of gene ids and names.
```
 eeisp_add_genename_from_geneid \
     Sample_CDI_score_data_thre0.5.txt \
     Sample_CDI_score_data_thre0.5.addgenename.txt \
     geneidlist.txt
 eeisp_add_genename_from_geneid \
     Sample_EEI_score_data_thre0.5.txt \
     Sample_EEI_score_data_thre0.5.addgenename.txt \
     geneidlist.txt
```
The output files include gene names.

```
   $ head Sample_CDI_score_data_thre0.5.addgenename.txt
   2       7       ESG000003       ESG000008       OR4F5   FO538757.3      0.96384320244841
   0       1       ESG000001       ESG000002       RP11-34P13.3    FAM138A 0.6852891560232545
   0       6       ESG000001       ESG000007       RP11-34P13.3    RP11-34P13.9    0.6852891560232545
   7       8       ESG000008       ESG000009       FO538757.3      FO538757.2      0.6852891560232545
   3       9       ESG000004       ESG000010       RP11-34P13.7    AP006222.2      0.6469554204484568
   4       6       ESG100005       ESG000007       RP11-34P13.8    RP11-34P13.9    0.5258703930217091
```

## Reference
Nakajima N., Hayashi T., Fujiki K., Shirahige K., Akiyama T., Akutsu T. and Nakato R., [Codependency and mutual exclusivity for gene community detection from sparse single-cell transcriptome data](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkab601/6324613), *Nucleic Acids Research*, 2021.
