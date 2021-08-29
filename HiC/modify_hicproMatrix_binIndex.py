#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan

'''
%prog <iced_matrix>

The bin index in HiC-Pro iced_matirx should start from 1 instead of 0
'''

import re
import sys
import pandas as pd


def modify_icedmatrix(matrix):
    prefix = re.findall('([0-9a-zA-Z\_\.\-]+)\.matrix',matrix_file)[0]
    df = pd.read_table(matrix,sep='\t',names=['bin1','bin2','score'])
    df['newbin1'] = df['bin1'].astype(int)+1
    df['newbin2'] = df['bin2'].astype(int)+1
    newdf = df.loc[:,['newbin1','newbin2','score']]
    newdf.to_csv('{}.modified.matrix'.format(prefix),sep='\t',header=False,index=False)

if __name__ == '__main__':
    from optparse import OptionParser

    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args) == 1:
        matrix_file = args[0]
        modify_icedmatrix(matrix_file)
    else:
        sys.exit(p.print_help())
