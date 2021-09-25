#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
 EEISP computes the score of CDI and EEI for 
 positive and negative samples from glioblastoma scRNA-seq data
'''

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
  # Count the number of genes and cells.
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
  data_file0 = str(filename) + '_number_nonzero_exp_gene.txt'
  #data_file0 = 'Stem_P10_number_nonzero_exp_gene.txt'
  fout0 = open(data_file0, "w")
  for d in range(0, len(count_exp)):
      fout0.writelines(str(d) +"\t"+ str(nzero_index[d])+"\t"+ str(count_exp[d])+"\n")
  fout0.close()
  print ("-----------------------------------------------")  


def generate_EEI(A, tarray, farray, filename):
  Allgene = A.shape[0]
  Allcell = A.shape[1]
  Count_excl = np.zeros((Allgene, Allgene), dtype=np.int64)
  is_nonzeroMat = A > 0             
  is_reverseMat = np.logical_not(A)
  for i in range(Allgene):
     Count_excl[i] = np.sum(np.logical_and(is_nonzeroMat[i], is_reverseMat), axis=1)
  print ("Count the number of cells that two genes are expressed exclusively----")
  data_file0 = str(filename) + '_data_exclusive.txt' 
  #np.savetxt(data_file0, Count_excl, delimiter="\t")

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

  total_array = calc_EEI_matrix( Allgene, Allcell, tarray, farray, Count_excl, Prob_excl )

  return total_array


def calc_EEI_matrix( Allgene, Allcell, tarray, farray, Count_excl, Prob_excl ):
   EEI = np.zeros((Allgene, Allgene), dtype=np.float64)
   arrayb = []
   arrayc = []
   for i in range(0, Allgene):
     for j in range(i+1, Allgene):
       for a1 in range(0, len(tarray)):
         if( (i == tarray[a1][0] and j == tarray[a1][1]) or (j == tarray[a1][0] and i == tarray[a1][1]) ):
           score_t = calc_EEI( i, j, a1, Allcell, tarray, Count_excl, Prob_excl)
           EEI[i][j] = EEI[j][i] = score_t           
           print ("EEI(%d,%d)=%.3F, EEI(%d,%d)=%.3F" % (i,j,EEI[i][j],j,i,EEI[j][i]))
           print ("-----------------------------------------------")
           arrayt = store_EEI( i, j, EEI[i][j], 1 )
           arrayb.append(arrayt)

       for a2 in range(0, len(farray)):
         if( (i == farray[a2][0] and j == farray[a2][1]) or (j == farray[a2][0] and i == farray[a2][1]) ):
           score_f = calc_EEI( i, j, a2, Allcell, farray, Count_excl, Prob_excl)
           EEI[i][j] = EEI[j][i] = score_f  
           print ("EEI(%d,%d)=%.3F, EEI(%d,%d)=%.3F" % (i,j,EEI[i][j],j,i,EEI[j][i]))
           print ("-----------------------------------------------")
           arrayf = store_EEI( i, j, EEI[i][j], 0 )
           arrayc.append(arrayf)

   df1 = pd.DataFrame(arrayb)
   df2 = pd.DataFrame(arrayc)
   df_array = pd.concat([df1, df2])

   return df_array


val = 0.0
def calc_EEI( num1, num2, p, ncell, array, Count_excl, Prob_excl):
   # Compute EEI for gene pairs.
   global val
   x1 = Count_excl[num2][num1]
   p1 = Prob_excl[num1][num2]
   prob1 = binom.sf(x1-1, ncell, p1)
   x2 = Count_excl[num1][num2]
   p2 = Prob_excl[num2][num1]
   prob2 = binom.sf(x2-1, ncell, p2)
   if( prob1 <= 0 or prob2 <= 0 ): 
      val = -10000000.0 
   else:
      val = ((-(math.log10(prob1))) + (-(math.log10(prob2)))) / 2 

   return val


def store_EEI( num1, num2, num_eei, cnum ): 
   gene1 = num1
   gene2 = num2
   eei = num_eei
   classify = cnum
   gene = []
   gene.append(gene1)   
   gene.append(gene2) 
   gene.append(eei)   
   gene.append(classify)
  
   return gene


def write_EEI( df, filename ):

   data_file0 = str(filename) + '_True_False_EEI_score_data.txt' 
   df.to_csv(data_file0, sep="\t", header=False, index=False)

              
def main():

  parser = argparse.ArgumentParser()
  parser.add_argument("matrix", help="Input matrix", type=str)
  parser.add_argument("filename", help="File name", type=str)
  parser.add_argument("positive", help="Positive samples", type=str)
  parser.add_argument("negative", help="Negative samples", type=str)
  args = parser.parse_args()
  print(args)
  
  startt = time.time()

  input_data = pd.read_csv(args.matrix, index_col=0, sep="\t")
  true_data = pd.read_table(args.positive, header=None)
  tarray = np.array(true_data)
  false_data = pd.read_table(args.negative, header=None)
  farray = np.array(false_data)
  
  A = count_genes(input_data, args.filename)
  
  array = generate_EEI(A, tarray, farray, args.filename)
 
  write_EEI(array, args.filename)  
 
  print("Finish to compute EEI score!")
  elapsed_time = time.time() - startt
  print ("Elapsed_time:{0}".format(elapsed_time) + "[sec]")
  print("*************************************************************")


if __name__ == '__main__':
    main()

