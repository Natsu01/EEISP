#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
EEISP computes the score of CDI and EEI from scRNA-seq data
'''

import numpy as np
import pandas as pd
import math
import sys
import time
from math import log
from scipy.stats import binom
import argparse
import multiprocessing as mp
from multiprocessing import Pool

def argwrapper(args):
    return args[0](*args[1:])

def calcCDI_eachrow(i, ProbArray, CountArray, ngene, ncell):
    array = np.zeros(ngene)
    for j in range(i + 1, ngene):
        x = CountArray[j]
        p = ProbArray[j]
        # 生存関数 binom.sfを使う
         # xが入っていないのでx-1に気をつけること
        prob = binom.sf(x - 1, ncell, p)
        if prob <= 0:
            val = -10000000.0
        else:
            val = -(math.log10(prob))
        array[j] = val

    del ProbArray
    del CountArray
    return array

'''
各遺伝子ペアで(gi,gj)=(0,1)となるサンプル数がCount_excl[j][i]回以上になるp値(3prob1)と排他性スコアの計算
各遺伝子ペアで(gi,gj)=(1,0)となるサンプル数がCount_excl[i][j]回以上になるp値(prob2)と排他性スコアの計算
便宜的に2つのp値の平均を計算する
'''
def calcEEI_eachrow(i, ProbArray_i, CountArray_i, ProbArray_T_i, CountArray_T_i, ngene, ncell):
    array = np.zeros(ngene)
    for j in range(i + 1, ngene):
        x1 = CountArray_T_i[j]
        p1 = ProbArray_i[j]
        prob1 = binom.sf(x1 - 1, ncell, p1)
        x2 = CountArray_i[j]
        p2 = ProbArray_T_i[j]
        prob2 = binom.sf(x2 - 1, ncell, p2)
        if prob1 <= 0 or prob2 <= 0:
            val = -10000000.0
        else:
            val = ((-(math.log10(prob1))) + (-(math.log10(prob2)))) / 2
        array[j] = val
    return array


def genMatrix_MultiProcess(Prob_joint, Count_joint, MatType, ngene, ncell, *, ncore=4):
    context = mp.get_context('spawn')

    p = Pool(ncore)
    func_args = []

    for i in range(0, ngene):
        if MatType == "CDI":
            func_args.append((calcCDI_eachrow, i, Prob_joint[i], Count_joint[i], ngene, ncell))
        elif MatType == "EEI":
            func_args.append((calcEEI_eachrow, i, Prob_joint[i], Count_joint[i], Prob_joint.T[i], Count_joint.T[i], ngene, ncell))
        else:
            print("Error: illegal MatType for genMatrix_MultiProcess.")
            sys.exit()

    results = p.map(argwrapper, func_args)
    p.close()

    Matrix = np.array(results)
    Matrix = Matrix + Matrix.T - np.diag(np.diag(Matrix))

    return Matrix


def generate_CDImatrix(A, args):
    ngene = A.shape[0]
    ncell = A.shape[1]
    is_nonzeroMat = A > 0
    p_nonzero = np.sum(is_nonzeroMat, axis=1) / ncell

    if(args.gpu):
        print("using GPU for CDI calculation.")
        import cupy as cp
        p_nonzero = cp.asarray(p_nonzero)
        is_nonzeroMat = cp.asarray(is_nonzeroMat)
        Prob_joint = p_nonzero * p_nonzero[:, cp.newaxis]
        Count_joint = cp.zeros((ngene, ngene), dtype=cp.int64)
        for i in range(ngene):
            Count_joint[i] = cp.sum(is_nonzeroMat[i] * is_nonzeroMat, axis=1)

        Prob_joint = cp.asnumpy(Prob_joint)
        Count_joint = cp.asnumpy(Count_joint)
    else:
        print("using CPU for CDI calculation.")
        Prob_joint = p_nonzero * p_nonzero[:, np.newaxis]
        Count_joint = []
        for row in is_nonzeroMat:
            Count_joint.extend(np.sum(row * is_nonzeroMat, axis=1))
        Count_joint = np.array(Count_joint).reshape(ngene, ngene)

#    np.savetxt(args.output + "_number_joint_nonzero_gene_thre10.0.txt", Count_joint, delimiter="\t")

    CDI = genMatrix_MultiProcess(Prob_joint, Count_joint, "CDI", ngene, ncell, ncore=args.threads)

    return CDI

def generate_EEImatrix(A, args):
    ngene = A.shape[0]
    ncell = A.shape[1]
    is_nonzeroMat = A > 0
    p_nonzero = np.sum(is_nonzeroMat, axis=1) / ncell
    p_zero = np.sum(A == 0, axis=1) / ncell
#    np.savetxt(args.output + "_prob_nonzero.txt", p_nonzero, delimiter="\t")
#    np.savetxt(args.output + "_prob_zero.txt", p_zero, delimiter="\t")

    if(args.gpu):
        print("using GPU for EEI calculation.")
        import cupy as cp
        p_nonzero = cp.asarray(p_nonzero)
        p_zero = cp.asarray(p_zero)
        is_nonzeroMat = cp.asarray(is_nonzeroMat)
        notA = cp.asarray(np.logical_not(A))

        Prob_joint = p_nonzero * p_zero[:, np.newaxis]
        Count_excl = cp.zeros((ngene, ngene), dtype=np.int64)
        for i in range(ngene):
            Count_excl[i] = cp.sum(cp.logical_and(is_nonzeroMat[i], notA), axis=1)

        Prob_joint = cp.asnumpy(Prob_joint)
        Count_excl = cp.asnumpy(Count_excl)
    else:
        print("using CPU for EEI calculation.")
        # 各遺伝子ペアで排他的発現をする確率の初期化と確率計算
        Prob_joint = p_nonzero * p_zero[:, np.newaxis]
        # 1回の行列演算で排他的発現(g1,g2)=(1,0) と(g2,g1)=(1,0)=(g1,g2)=(0,1)を持つサンプル数をカウント
        notA = np.logical_not(A)
        Count_excl = []
        for row in is_nonzeroMat:
            # Aと(NOT A)の転置行列をかけた行列の、各要素のサンプル数をカウント
            Count_excl.extend(np.sum(np.logical_and(row, notA), axis=1))
        Count_excl = np.array(Count_excl).reshape(ngene, ngene)

#    np.savetxt(args.output + "_data_exclusive.txt", Count_excl, delimiter="\t")

    EEI = genMatrix_MultiProcess(Prob_joint, Count_excl, "EEI", ngene, ncell, ncore=args.threads)
    return EEI


def calc_degree(Matrix, thre, ngene, filename, output, genenames):
    df = pd.DataFrame(Matrix)
    degree = np.sum(df > thre).tolist()
    df = df[df > thre]
    df = df.stack().reset_index()
    df.columns = ["i", "j", "val"]
    df = df[df["i"] < df["j"]]
    df["gene_i"] = genenames[df["i"]]
    df["gene_j"] = genenames[df["j"]]
    df = df.reindex(columns=["i", "j", "gene_i", "gene_j", "val"])
    df = df.sort_values(["val", "i"], ascending=[False, True])

    data_file = output + "_" + filename + ".txt"
    print ("output degree data in " + data_file)
    print ("number of gene pairs over threshold (>" + str(thre) + "): " + str(df.shape[0]))
    df.to_csv(data_file, sep="\t", header=False, index=False)
    return degree


def calc_degree_dist(degree, filename, args):
    max_value = max(degree)
    min_value = min(degree)
    value_range = max_value - min_value
    print("max degree:%.3F min degree:%.3F value_width=%.3F" % (max_value, min_value, value_range))

    freq = []
    for a in range(min_value + 1, max_value + 1):
        fnum = degree.count(a)
        if fnum > 0:
            freq.append([a, fnum])

    df = pd.DataFrame(freq, columns=["Degree", "The number of genes"])

    log_df = np.log(df)
    log_df = log_df.rename(
        columns={
            "Degree": "Log_Degree",
            "The number of genes": "Log_The number of genes",
        }
    )
    merge = pd.concat([log_df, df], axis=1)
    merge.to_csv(args.output + "_" + filename + "_degree_distribution.tsv", sep="\t")


def get_nonzero_matrix(input_data, output):
    A = np.array(input_data)
    A = A[np.any(A > 0, axis=1)]

    ncell_exp = np.sum(input_data > 0, axis=1)
    df = pd.DataFrame(ncell_exp[ncell_exp>0])
    df["id"] = range(df.shape[0])
    df["genename"] = df.index
    df["count"] = df[0]
    df[["id","genename","count"]].to_csv(output + '_number_nonzero_exp.txt', sep="\t", index=False, header=False)

    genenames = df.index
    return A, genenames


def main():
    parser = argparse.ArgumentParser(prog='eeisp')
    parser.add_argument("matrix", help="Input matrix", type=str)
    parser.add_argument("output", help="Output prefix", type=str)
    parser.add_argument("--threCDI", help="Threshold for CDI (default: 20.0)", type=float, default=20)
    parser.add_argument("--threEEI", help="Threshold for EEI (default: 10.0)", type=float, default=10)
    parser.add_argument("--tsv", help="Specify when the input file is tab-delimited (.tsv)", action="store_true")
    parser.add_argument("--gpu", help="GPU mode", action="store_true")
    parser.add_argument("-p", "--threads", help="number of threads (default: 2)", type=int, default=2)
    parser.add_argument("-v", "--version", action='version', version='%(prog)s 0.2.0')

    args = parser.parse_args()
    print(args)

    startt = time.time()

    if (args.tsv):
        input_data = pd.read_csv(args.matrix, index_col=0, sep="\t")
    else:
        input_data = pd.read_csv(args.matrix, index_col=0, sep=",")

    print("number of cells: ", input_data.shape[1])
    print("number of genes: ", input_data.shape[0])

    A, genenames = get_nonzero_matrix(input_data, args.output)
    ngene = A.shape[0]
    ncell = A.shape[1]
    print("number of nonzero genes: ", ngene)
    print ("-----------------------------------------------")

    CDI = generate_CDImatrix(A, args)
    EEI = generate_EEImatrix(A, args)
    DEGREE_CDI = calc_degree(CDI, args.threCDI, ngene, "CDI_score_data_thre" + str(args.threCDI), args.output, genenames)
    DEGREE_EEI = calc_degree(EEI, args.threEEI, ngene, "EEI_score_data_thre" + str(args.threEEI), args.output, genenames)

    print("Finish Co-dependency Network!")
    elapsed_time = time.time() - startt
    print("Elapsed_time:{0}".format(elapsed_time) + "[sec]")
    print("*************************************************************")

    # CDI, EEIの次数分布
    calc_degree_dist(DEGREE_CDI, "CDI", args)
    calc_degree_dist(DEGREE_EEI, "EEI", args)

    print("Finish to write the CDI and EEI degreee distribution!")


if __name__ == "__main__":
    main()
