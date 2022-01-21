#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2022/1/16

'''
%prog <GeneCount file> <OGvenn.csv>

Counting the number of all genes according to OGids
OGnumberfile.txt : according to jvenn website manual statistics to get

>>> python %prog Orthogroups.GeneCount.tsv OGvenn.csv

__author__ = Haoran Pan
__mail__ = panpyhr@gmail.com
__date__ = 20220121
__version__ = v1.1
'''

import sys
import pandas as pd
from rich import print


def read_Genecount(file):
    df = pd.read_table(file, sep='\t')
    return df


def read_OGcsv(OGcsv):
    OGdf = pd.read_csv(OGcsv)
    return OGdf


def stat_OGnum_genes(GeneCdf,OGdf):
    OGD = {}
    col_name = [column for column in OGdf]
    for col in col_name:
        qOG = [i for i in OGdf[col].tolist() if pd.isnull(i) == False]
        qdf = GeneCdf[GeneCdf['Orthogroup'].isin(qOG)]
        # print(qdf)
        #print(qdf.shape[0],len(qOG))
        if qdf.shape[0] == len(qOG):
            Gtotal = qdf['Total'].sum()
            # print(Gtotal)
            OGD[col] = [len(qOG),Gtotal]
        else:
            print('[bold magenta]{} rows != len(OGids) !!! Please check it ![/bold magenta]'.format(col))
    return OGD


def output_OGD(OGD):
    for key in OGD.keys():
        print('[bold green]{}[/bold green]\t[bold yellow]{}[/bold yellow]\t[bold red]{}[/bold red]'.format(key,OGD[key][0],OGD[key][1]))


if __name__ == '__main__':
    from optparse import OptionParser
    from rich.traceback import install
    install()
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args) == 2:
        count_file = args[0]
        OGcsv_file = args[1]
        df = read_Genecount(count_file)
        OGdf = read_OGcsv(OGcsv_file)
        OGD = stat_OGnum_genes(df,OGdf)
        output_OGD(OGD)
    else:
        sys.exit(p.print_help())
