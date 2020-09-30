#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# @Time : 2020/8/19 10:48


import pandas as pd
import argparse

"""
此脚本用于在表达谱中提取目的基因的表达数据
"""

def gain_gene_expression_in_matrix(expression_matrix,geneID_file,output_file):
    if 'xlsx' or 'xls' in expression_matrix:
        df = pd.read_excel(expression_matrix)
    elif 'csv' in expression_matrix:
        df = pd.read_csv(expression_matrix)
    else:
        df = pd.read_table(expression_matrix,sep='\t',header=0)
    with open(geneID_file) as f:
        geneID_lsit = [line.strip() for line in f]
    concat_list = []
    for index,row in df.iterrows():
        if row[0] in geneID_lsit:
            concat_list.append(df.loc[[index],:])
    result = pd.concat(concat_list)
    result.to_excel(output_file,header=True,index=False)

def main(args):
    expression_matirx = args.matrix
    geneID_file = args.geneid
    output_file = args.output
    gain_gene_expression_in_matrix(expression_matirx,geneID_file,output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='Extract_gene_expression_in_matrix.py',
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Extract the expression data of the target gene from the expression profile''')
    parser.add_argument('-m', '--matrix', required=True, help='Input the expression profile file in the RNA-Seq results(.xlsx)')
    parser.add_argument('-g', '--geneid', required=True, help='Input the target gene file')
    parser.add_argument('-o', '--output', required=True, help='Give the output file a name(.xlsx)')
    args = parser.parse_args()
    main(args)
