#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import csv
import math
import argparse


def main():

  parser = argparse.ArgumentParser()
  parser.add_argument("filename", help="filename", type=str)
  parser.add_argument("file0", help="gene_file", type=str)
  parser.add_argument("file1", help="cdi_score_file", type=str)
  parser.add_argument("file2", help="eei_score_file", type=str)
  parser.add_argument("file3", help="expression_file", type=str)
  parser.add_argument("thre", help="threshold", type=str)

  args = parser.parse_args()
  print(args)

  input_data0 = pd.read_csv(args.file0, dtype='str', header=None, sep="\t").values.tolist()
  size0 = len(input_data0) 

  input_data1 = pd.read_csv(args.file1, dtype='str', header=None, sep="\t").values.tolist()
  size1 = len(input_data1)

  input_data2 = pd.read_csv(args.file2, dtype='str',  header=None, sep="\t").values.tolist()
  size2 = len(input_data2)

  input_data3 = pd.read_csv(args.file3, dtype='str',  header=None, sep="\t").values.tolist()
  size3 = len(input_data3)


  node0 = []
  for i in range(0, size3):
    genenum = input_data3[i][0]
    geneID = input_data3[i][1]
    geneEXP = input_data3[i][2]
    rdata = []
    for j in range(0, size0):
      if( genenum == input_data0[j][0] ):	
       name = input_data0[j][2]
       rdata.append( genenum )
       rdata.append( geneID )
       rdata.append( name )
       rdata.append( geneEXP )
       
    if( len(rdata) == 4 ):
      node0.append(rdata)


  node1 = []
  for i in range(0, size1):
    g1 = input_data1[i][0]
    g2 = input_data1[i][1] 
    score = input_data1[i][2] 
    gene = []
    for j in range(0, len(node0)):
       if( g1 == node0[j][0] ):
         gene_name1 = node0[j][2]
         gene.append( gene_name1 )
         for k in range(0, len(node0)):
            if( g2 == node0[k][0] ):
               gene_name2 = node0[k][2]
               gene.append( gene_name2 )
               gene.append( score )
    
    if( len(gene) == 3 ):
       node1.append(gene)

  node2 = []
  for i in range(0, size2):
    g1 = input_data2[i][0]
    g2 = input_data2[i][1]
    score = input_data2[i][2]
    gene = []
    for j in range(0, len(node0)):
       if( g1 == node0[j][0] ):
         gene_name1 = node0[j][2]
         gene.append( gene_name1 )
         for k in range(0, len(node0)):
            if( g2 == node0[k][0] ):
               gene_name2 = node0[k][2]
               gene.append( gene_name2 )
               gene.append( score )

    if( len(gene) == 3 ):
       node2.append(gene)


  flag1 = np.zeros((size1, size1), dtype=np.int64)
  flag2 = np.zeros((size2, size2), dtype=np.int64)

  data_file_cdi = str(args.filename) + '_CDI_convert_data_thre' + str(args.thre) + '.txt'
  data_file_eei = str(args.filename) + '_EEI_convert_data_thre' + str(args.thre) + '.txt'
  
  fout1 = open(data_file_cdi, "w")
  count1 = 0
  for i in range(0, len(node1)):
    n1 = int(input_data1[i][0])
    n2 = int(input_data1[i][1])
    if( flag1[n1][n2] == 0 or flag1[n2][n1] == 0 ):
       fout1.writelines(str(node1[i][0]) +"\t"+ str(node1[i][1])+"\t"+ str(node1[i][2])+"\n")
       flag1[n1][n2] = flag1[n2][n1] = 1
       count1 = count1+1
  fout1.close()

  count2 = 0
  fout2 = open(data_file_eei, "w")
  count = 0
  for i in range(0, len(node2)):
    n1 = int(input_data2[i][0])
    n2 = int(input_data2[i][1])
    if( flag2[n1][n2] == 0 or flag2[n2][n1] == 0 ):
       fout2.writelines(str(node2[i][0]) +"\t"+ str(node2[i][1])+"\t"+ str(node2[i][2])+"\n")
       flag2[n1][n2] = flag2[n2][n1] = 1
       count2 = count2+1
  fout2.close()


  print("*************************************************************");
  print("Finish to convert to gene names!");
  print("The number of CDI gene pairs : %d"  %(count1));
  print("The number of EEI gene pairs : %d"  %(count2));
  print("*************************************************************"); 

       
if __name__ == '__main__':
    main()

   
