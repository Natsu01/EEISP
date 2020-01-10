import numpy as np
import pandas as pd
import csv
import argparse



def main():

  parser = argparse.ArgumentParser()
  parser.add_argument("filename", help="filename", type=str)
  parser.add_argument("file1", help="score_file", type=str)
  parser.add_argument("file2", help="expression_file", type=str)
  parser.add_argument("score", help="score", type=str)
  parser.add_argument("thre", help="threshold", type=str)

  args = parser.parse_args()
  print(args)
  
  node = []
  input_data1 = pd.read_csv(args.file1, header=None, sep="\t").values.tolist()
  size1 = len(input_data1)
  #print(input_data1.index) 

  input_data2 = pd.read_csv(args.file2, header=None, sep="\t").values.tolist()
  size2 = len(input_data2)
  #print(input_data2.index) 
  
  for i in range(0, size1):
    g1 = input_data1[i][0]
    g2 = input_data1[i][1] 
    cdi = input_data1[i][2] 
    gene = []
    for j in range(0, size2):
       if( g1 == input_data2[j][0] ):
         gene_name1 = input_data2[j][1]
         gene.append( gene_name1 )
         for k in range(0, size2):
            if( g2 == input_data2[k][0] ):
               gene_name2 = input_data2[k][1]
               gene.append( gene_name2 )
               gene.append( cdi )
    
    if( len(gene) == 3 ):
       node.append(gene)


  flag = np.zeros((size1, size1), dtype=np.int64)
  num = int(args.number)
  if( num == 0 ):
    data_file0 = str(args.filename) + '_CDI_convert_data_thre' + str(args.thre) + '.txt'
  else:
    data_file0 = str(args.filename) + '_EEI_convert_data_thre' + str(args.thre) + '.txt'
  
  fout0 = open(data_file0, "w")
  count = 0
  for i in range(0, len(node)):
    n1 = int(input_data1[i][0])
    n2 = int(input_data1[i][1])
    if( flag[n1][n2] == 0 or flag[n2][n1] == 0 ):
       fout0.writelines(str(node[i][0]) +"\t"+ str(node[i][1])+"\t"+ str(node[i][2])+"\n")
       flag[n1][n2] = flag[n2][n1] = 1
       count = count+1
  fout0.close()   


  print("*************************************************************");
  print("Finish to　write to the file!");
  print("The number of all gene pairs : %d"  %(count));
  print("*************************************************************"); 

       
if __name__ == '__main__':
    main()

   
