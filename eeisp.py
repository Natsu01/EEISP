# EEISP computes the score of CDI and EEI from scRNA-seq data


import numpy as np
import pandas as pd
import csv
import decimal
import math
import time
from math import log
from math import ceil
import sys
import argparse
from scipy.stats import binom
from operator import itemgetter


def count_genes(input_data, filename):
  # Count the number of genes and cells
  A = np.array(input_data)
  Allgene = A.shape[0]
  Allcell = A.shape[1]
  print ("Matrix size of A:", A.shape)
  print ("Number of all genes:", Allgene) 
  print ("Number of all cells:", Allcell)
 
  zero = np.all(A == 0, axis=1)                
  notzero = input_data[np.logical_not(zero)]            
  nzero_index = notzero.index
  print(nzero_index)
  
  A = A[np.logical_not(np.all(A == 0, axis=1))]
  print(A.shape)
  Allgene = A.shape[0]
  print ("Number of nonzero genes:", Allgene)
  print ("Number of all cells:", Allcell)

  count_cells(A, nzero_index, filename)

  return A


def count_cells(A, nzero_index, filename):
  # Count the number of cells that have nonzero expression for each gene 
  # and store the gene number, EnsembleID(gene name) and the number of cells 
  # that have nonzero expression into array 'count_exp'. 
  count_exp = np.sum(A > 0, axis=1)
  data_file0 = str(filename) + '_number_nonzero_exp.txt'
  fout0 = open(data_file0, "w")
  for d in range(0, len(count_exp)):
      fout0.writelines(str(d) +"\t"+ str(nzero_index[d])+"\t"+ str(count_exp[d])+"\n")
  fout0.close()
  print ("-----------------------------------------------")  


def calc_degree_CDI(A, threCDI, filename ):
  Allgene = A.shape[0]
  Allcell = A.shape[1]
  # Count the number of cells that have nonzero expression for each gene 
  # and compute the probability that a gene have nonzero expression.
  p_nonzero = np.sum(A > 0, axis=1) / Allcell

  # Initialize the matrix of 'Prob_joint' and compute the joint probability 
  # that two genes gi and gj have jointly nonzero expression (gi=1,gj=1). 
  Prob_joint = p_nonzero * p_nonzero[:, np.newaxis]

  # Count the number of cells that two genes jointly have nonzero expression.  
  Count_joint = np.zeros((Allgene, Allgene), dtype=np.int64)
  is_nonzeroMat = A > 0
  for i in range(Allgene):
     Count_joint[i] = np.sum(is_nonzeroMat[i] * is_nonzeroMat, axis=1)
  data_file0 = str(filename) + '_joint_nonzero_gene_thre' + str(threCDI) + '.txt'
  #np.savetxt(data_file0, Count_joint, delimiter='\t')
  print ("Count the number of cells that two genes are jointly expressed----")  

  return Count_joint, Prob_joint



def calc_degree_EEI(A, threEEI, filename ):
  Allgene = A.shape[0]
  Allcell = A.shape[1]
  # Count the number of cells that exhibit the mutually exclusive expression, 
  # (gi,gj)=(1,0) and (gi,gj)=(0,1). 
  Count_excl = np.zeros((Allgene, Allgene), dtype=np.int64)
  is_nonzeroMat = A > 0
  is_reverseMat = np.logical_not(A)
  for i in range(Allgene):
     Count_excl[i] = np.sum(np.logical_and(is_nonzeroMat[i], is_reverseMat), axis=1)
  data_file0 = str(filename) + '_data_exclusive.txt' 
  #np.savetxt(data_file0, Count_excl, delimiter='\t')
  print ("Count the number of cells that two genes are expressed exclusively----")

  # Count the number of cells that have nonzero expression and zero expression
  # for each gene and compute the probability of them. 
  p_nonzero = np.sum(A > 0, axis=1) / Allcell
  p_zero = np.sum(A == 0, axis=1) / Allcell
  data_file1 = str(filename) + '_prob_nonzero.txt' 
  data_file2 = str(filename) + '_prob_zero.txt'
  #np.savetxt(data_file1, p_nonzero, delimiter='\t')
  #np.savetxt(data_file2, p_zero, delimiter='\t')

  # Initialize a matrix of 'Prob_excl' and compute the probability that 
  # each gene pair exhibits the mutually exclusive expression. 
  Prob_excl = p_nonzero * p_zero[:, np.newaxis]
  data_file3 = str(filename) + '_data_prob_exlusive.txt'
  #np.savetxt(data_file3, Prob_excl, delimiter='\t')

  return Count_excl, Prob_excl



def calc_CDI(A, Count_joint, Prob_joint, Count_excl, Prob_excl, threCDI, threEEI, filename):
  # Compute the probability that the number of cells that two genes jointly have 
  # nonzero expression is more than 1_{ij} using a binomial distribution and the score of CDI. 
  Allgene = A.shape[0]
  Allcell = A.shape[1]
  CDI = np.zeros((Allgene, Allgene), dtype=np.float64)
  EEI = np.zeros((Allgene, Allgene), dtype=np.float64)
  DEGREE_CDI = []
  DEGREE_EEI = []
  cdineighbor = [] 
  eeineighbor = []
  CDI_neighbor = []
  EEI_neighbor = []
  degree_cdi = []
  degree_eei = []
  for i in range(0, Allgene):
     for j in range(i+1, Allgene):
       x = Count_joint[i][j]
       p = Prob_joint[i][j]
       n = Allcell
       prob = binom.sf(x-1, n, p)
       if( prob <= 0 ):
          CDI[i][j] = CDI[j][i] = -10000000.0 
       else:
          CDI[i][j] = CDI[j][i] = -(math.log10(binom.sf(x-1, n, p)))
       #print ("CDI(%d,%d)=%.3F, CDI(%d,%d)=%.3F" % (i,j,CDI[i][j],j,i,CDI[j][i]))
       #print ("-----------------------------------------------") 
  
       x1 = Count_excl[j][i]
       p1 = Prob_excl[i][j]
       n = Allcell
       prob1 = binom.sf(x1-1, n, p1)
       x2 = Count_excl[i][j]
       p2 = Prob_excl[j][i]
       prob2 = binom.sf(x2-1, n, p2)
       if( prob1 <= 0 or prob2 <= 0 ): 
          EEI[i][j] = EEI[j][i] = -10000000.0 
       else:
          score_eei = ((-(math.log10(prob1))) + (-(math.log10(prob2)))) / 2  
          EEI[i][j] = EEI[j][i] = score_eei  
       
       #print ("EEI(%d,%d)=%.3F, EEI(%d,%d)=%.3F" % (i,j,EEI[i][j],j,i,EEI[j][i]))
       #print ("-----------------------------------------------")

     array_cdi = []  
     array_eei = []
    
     for j in range(0, Allgene):
       if( CDI[i][j] != 0.0 or CDI[j][i] != 0.0 ):
          gene1 = i
          gene2 = j
          cdi = CDI[i][j]
          gene = []
          gene.append(gene1)  
          gene.append(gene2) 
          gene.append(cdi)   
          array_cdi.append(gene)
      
       if( EEI[i][j] != 0.0 or EEI[j][i] != 0.0 ):
          gene1 = i
          gene2 = j
          eei = EEI[i][j]
          gene = []
          gene.append(gene1)  
          gene.append(gene2) 
          gene.append(eei)    
          array_eei.append(gene)

     cdineighbor, cdineighbor_gene = sort_neighbor(i, array_cdi, 0, threCDI)
     CDI_neighbor.append(cdineighbor)
     degree_cdi = len(cdineighbor_gene)
     DEGREE_CDI.append(degree_cdi)
 
     eeineighbor, eeineighbor_gene = sort_neighbor(i, array_eei, 1, threEEI)
     EEI_neighbor.append(eeineighbor)
     degree_eei = len(eeineighbor_gene)
     DEGREE_EEI.append(degree_eei)

  write_file(CDI_neighbor, threCDI, 0, filename)
  write_file(EEI_neighbor, threEEI, 1, filename)

  return DEGREE_CDI, DEGREE_EEI



def sort_neighbor(i, array, num, threshold):
   if( num == 0 ):
     #print("Sort the co-expressed gene for %d-th node" % (i))
     array2 = sorted(array, key=itemgetter(2), reverse=True)
     print("The gene sets that CDI > %f for %d-th node" % (threshold,i))
   else:
     #print("Sort the mutually exclusive gene for %d-th node" % (i))
     array2 = sorted(array, key=itemgetter(2), reverse=True)
     print("The gene sets that EEI > %f for %d-th node" % (threshold,i))
  
   G = np.array(array2)
   print(G)
   neighbor = []
   neighbor_gene = []
   for a in range(0, G.shape[0]):
     ng = []
     if( G[a][2] > threshold ):  
         ng.append(i)  
         ng.append(int(G[a][1]))
         ng.append(float(G[a][2]))
         neighbor.append(ng)
         neighbor_gene.append(int(G[a][1]))

   return neighbor, neighbor_gene 



def write_file(neighbor, threshold, num, filename):
   if( num == 0 ):
      data_file0 = str(filename) + '_CDI_score_data_thre' + str(threshold) + '.txt' 
   else:
      data_file0 = str(filename) + '_EEI_score_data_thre' + str(threshold) + '.txt'

   fout0 = open(data_file0, "w")
   for i in range(0, len(neighbor)): 
      for d in range(0, len(neighbor[i])):
         fout0.writelines(str(neighbor[i][d][0]) +"\t"+ str(neighbor[i][d][1])+"\t"+ str(neighbor[i][d][2])+ "\n")
   fout0.close()




def calc_degree_dist(DEGREE, filename, num, threshold, dataname):
# Compute the degree distribution
  print(DEGREE)
  max_value = max(DEGREE)
  min_value = min(DEGREE)
  value_width = max_value - min_value
  
  if (num == 0 ):
     print ("max CDI degree:%.3F min CDI degree:%.3F value_width=%.3F" % (max_value, min_value, value_width))
  else: 
     print ("max EEI degree:%.3F min EEI degree:%.3F value_width=%.3F" % (max_value, min_value, value_width))

  freq = []
  for a in range(min_value+1, max_value+1):
     fnum = DEGREE.count(a)
     if (fnum > 0):
        freq.append([a, fnum])

  if (num == 0 ):
     df = pd.DataFrame(freq, columns=["CDI_degree", "The number of genes"])
     log_df = np.log(df)
     log_df = log_df.rename(columns={"CDI_degree": "Log_CDI_degree", "The number of genes": "Log_The number of genes"})
  else: 
     df = pd.DataFrame(freq, columns=["EEI_degree", "The number of genes"])
     log_df = np.log(df)
     log_df = log_df.rename(columns={"EEI_degree": "Log_EEI_degree", "The number of genes": "Log_The number of genes"})

  merge = pd.concat([log_df, df], axis=1)
  print(merge)

  if (num == 0): 
   merge.to_csv( str(dataname) + "_" + filename + "_thre" + str(threshold) + ".csv", sep='\t')
  else:
   merge.to_csv( str(dataname) + "_" + filename + "_thre" + str(threshold) + ".csv", sep='\t')


              
def main():

  parser = argparse.ArgumentParser()
  parser.add_argument("matrix", help="Input matrix", type=str)
  parser.add_argument("filename", help="File name", type=str)
  parser.add_argument("--threCDI", help="threshold for CDI (default: 20.0)", type=float, default=20)
  parser.add_argument("--threEEI", help="threshold for EEI (default: 10.0)", type=float, default=10)
  parser.add_argument("--version", help="estimateEEI version 1.0", type=float )

  args = parser.parse_args()
  print(args)
  
  startt = time.time()

  input_data = pd.read_csv(args.matrix, index_col=0, sep=",")
  print(input_data.index)  

  A = count_genes(input_data, args.filename)
  
  Count_joint, Prob_joint = calc_degree_CDI(A, args.threCDI, args.filename )
  Count_excl, Prob_excl = calc_degree_EEI(A, args.threEEI, args.filename )
  degree_cdi, degree_eei = calc_CDI(A, Count_joint, Prob_joint, Count_excl, Prob_excl, args.threCDI, args.threEEI, args.filename)
 
  print("--------------------------------------------------------------")
  print("Compute CDI and EEI degree distribution.")
  calc_degree_dist(degree_cdi, "CDI_degree_distribution", 0, args.threCDI, args.filename)
  calc_degree_dist(degree_eei, "EEI_degree_distribution", 1, args.threEEI, args.filename)  
  print("Finish to write the CDI and EEI degreee distribution!")
  print("--------------------------------------------------------------")
  print("Finish to compute the CDI and EEI scores!")
  elapsed_time = time.time() - startt
  print ("Elapsed_time:{0}".format(elapsed_time) + "[sec]")
  print("*************************************************************")



if __name__ == '__main__':
    main()

