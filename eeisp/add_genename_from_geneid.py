#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import sys
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input matrix", type=str)
    parser.add_argument("output", help="Output prefix", type=str)
    parser.add_argument("genelist", help="Gene list", type=str)
    parser.add_argument("--i_id", help="column number of gene id (default: 0)", type=int, default=0)
    parser.add_argument("--i_name", help="column number of gene name (default: 1)", type=int, default=1)

    args = parser.parse_args()

    i = args.i_id
    j = args.i_name
    genes = pd.read_csv(args.genelist, sep="\t", header=None)
    genes.index = genes[i]

    input_data = pd.read_csv(args.input, sep="\t", header=None)
    input_data.index = input_data[2]
    input_data["gene1"] = genes[j][input_data.index]
    input_data.index = input_data[3]
    input_data["gene2"] = genes[j][input_data.index]
    input_data.columns = ["i", "j", "geneid1","geneid2","val","genename1", "genename2"]
    input_data = input_data.reindex(columns=["i", "j", "geneid1","geneid2" ,"genename1", "genename2" ,"val"])

    input_data.to_csv(args.output, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
