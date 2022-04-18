#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan

'''
%prog <matirx file> <row num> <output.noctg.matrix>

>>> python %prog ZG_100000_iced.matrix 63228

Delete the contig level interaction signal from the Hi-C matrix file
row num : the index value corresponding to the last bin of the last chromosome
'''

import re
import sys
import pandas as pd

def remove_ctg_matrix(iced_matrix,rownum):
    prefix = re.findall('([0-9a-zA-Z\_\.\-]+)\.matrix', iced_matrix)[0]
    df = pd.read_table(iced_matrix, sep='\t', names=['bin1', 'bin2', 'score'])

    df = df[df['bin1'] <= int(rownum)]
    df = df[df['bin2'] <= int(rownum)]

    df.to_csv('{}.noctg.matrix'.format(prefix),sep='\t',header=False,index=False)



if __name__ == "__main__":
    from optparse import OptionParser
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args)==2:
        ks_file = args[0]
        remove_ctg_matrix(sys.argv[1],sys.argv[2])
    else:
        sys.exit(p.print_help())
