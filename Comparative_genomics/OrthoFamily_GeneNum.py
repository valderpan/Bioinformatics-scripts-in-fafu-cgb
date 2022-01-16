#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2022/1/16

'''
%prog <OGnumber file>

Counting the number of all genes according to OGids
OGnumberfile.txt : according to jvenn website manual statistics to get

>>> python %prog OGnumberfile.txt

__author__ = Haoran Pan
__mail__ = panpyhr@gmail.com
__date__ = 20220116
__version__ = v1.0
'''

import pandas as pd
import sys
from rich import print

def read_Genecount(file):
    df = pd.read_table(file, sep='\t')
    return df


def read_OGnumber(file):
    OGnumL = []
    with open(file) as f:
        lines = (line.strip() for line in f)
        for line in lines:
            if not line.startswith('OG'):
                continue
            else:
               OGnumL.append(line)
    return OGnumL


def stat_OGnum_genes(df,OGnumL):
    qdf = df[df['Orthogroup'].isin(OGnumL)]
    if qdf.shape[0] == len(OGnumL):
        Gtotal = qdf['Total'].sum()
        print(Gtotal)
    else:
        print('qdf rows != len(OGids) !!! Please check it !')


if __name__ == '__main__':
    from optparse import OptionParser
    from rich.traceback import install
    install()
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args) == 2:
        count_file = args[0]
        OGnum_file = args[1]
        df = read_Genecount(count_file)
        OGnumL = read_OGnumber(OGnum_file)
        stat_OGnum_genes(df,OGnumL)
    else:
        sys.exit(p.print_help())