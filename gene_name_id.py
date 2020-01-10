import numpy as np
import pandas as pd
import csv
import math
import argparse


def main():

  parser = argparse.ArgumentParser()
  parser.add_argument("filename", help="filename", type=str)
  parser.add_argument("file0", help="gene_file", type=str)
  parser.add_argument("file1", help="score_file", type=str)
  parser.add_argument("file2", help="expression_file", type=str)
  parser.add_argument("score", help="score", type=str)
  parser.add_argument("thre", help="threshold", type=str)

  args = parser.parse_args()
  print(args)

  input_data0 = pd.read_csv(args.file0, dtype='str', header=None, sep="\t").values.tolist()
  size0 = len(input_data0)
  #print(input_data0.index) 

  input_data1 = pd.read_csv(args.file1, dtype='str', header=None, sep="\t").values.tolist()
  size1 = len(input_data1)
  #print(input_data1) 

  input_data2 = pd.read_csv(args.file2, dtype='str',  header=None, sep="\t").values.tolist()
  size2 = len(input_data2)
  #print(input_data2.index) 


  node0 = []
  for i in range(0, size2):
    genenum = input_data2[i][0]
    geneID = input_data2[i][1]
    geneEXP = input_data2[i][2]
    rdata = []
    for j in range(0, size0):
      if( geneID == input_data0[j][1] ):	
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
    cdi = input_data1[i][2] 
    gene = []
    for j in range(0, len(node0)):
       if( g1 == node0[j][0] ):
         gene_name1 = node0[j][2]
         gene.append( gene_name1 )
         for k in range(0, len(node0)):
            if( g2 == node0[k][0] ):
               gene_name2 = node0[k][2]
               gene.append( gene_name2 )
               gene.append( cdi )
    
    if( len(gene) == 3 ):
       node1.append(gene)


  flag = np.zeros((size1, size1), dtype=np.int64)
  num = int(args.score)
  if( num == 0 ):
    data_file0 = str(args.filename) + '_CDI_convert_data_thre' + str(args.thre) + '.txt'
  else:
    data_file0 = str(args.filename) + '_EEI_convert_data_thre' + str(args.thre) + '.txt'
  
  fout0 = open(data_file0, "w")
  count = 0
  for i in range(0, len(node1)):
    n1 = int(input_data1[i][0])
    n2 = int(input_data1[i][1])
    if( flag[n1][n2] == 0 or flag[n2][n1] == 0 ):
       fout0.writelines(str(node1[i][0]) +"\t"+ str(node1[i][1])+"\t"+ str(node1[i][2])+"\n")
       flag[n1][n2] = flag[n2][n1] = 1
       count = count+1
  fout0.close()   


  print("*************************************************************");
  print("Finish to　write to the file!");
  print("The number of all gene pairs : %d"  %(count));
  print("*************************************************************"); 

       
if __name__ == '__main__':
    main()

   
