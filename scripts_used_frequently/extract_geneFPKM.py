#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/8/29


import re
import sys
import pandas as pd
import richlog as rl

__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20210829'
__version__ = 'v1.0'


def read_gene(genefile):
    G_D = {}
    for gene in genefile:
        G_L = []
        with open(gene,'r') as f:
            lines = (line.strip() for line in f)
            for line in lines:
                G_L.append(line)
        gene_p = re.findall('([0-9A-Za-z\.\_\-]+)\.txt',gene)[0]
        G_D[gene_p] = G_L
    return G_D


def read_fpkm(file):
    files = {}
    for f in file:
        fpkm_p = re.findall('([0-9A-Za-z\.\_\-]+)\.xlsx',f)[0]
        df = pd.read_excel(f)
        files[fpkm_p] = df
    return files


def ExtractOut_fpkm(geneD,files):
    for fpkm in files.keys():
        rl.info_out('Read the expression file : {}.xlsx'.format(fpkm))
        for gene in geneD.keys():
            rl.info_out('Read geneID file : {}.txt'.format(gene))
            qdf = files[fpkm][files[fpkm]['gene_id'].isin(geneD[gene])]
            rl.info_out('Output result file : {}-{}.xlsx'.format(fpkm,gene))
            qdf.to_excel('{}-{}.xlsx'.format(fpkm,gene),header=True,index=False)



def main(args):
    genefiles = args.gene
    fpkmfiles = args.fpkm
    geneD = read_gene(genefiles)
    files = read_fpkm(fpkmfiles)
    ExtractOut_fpkm(geneD,files)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Extraction of fpkm for a specific gene set\nNote1: The column name of the gene name column must be "gene_id"\nNote2: The suffix of geneID file must be [.txt]; the suffix of fpkm file must be [.xlsx]''',
        usage="python {} -g genefiles -f fpkmfiles ".format(sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,
                                                                            __mail__, __date__, __version__))
    parser.add_argument('-g', '--gene', required=True, nargs='+',
                        help='Input the geneID file to be extracted, you can input more than one')
    parser.add_argument('-f', '--fpkm', required=True, nargs='+',
                        help='Input the fpkm file, you can input more than one')
    args = parser.parse_args()
    main(args)
