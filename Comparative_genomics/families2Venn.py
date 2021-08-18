#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/4/18

'''
python %prog <GeneCount file>

Gene family clustering based on orthofinder results

__author__ = Haoran Pan
__mail__ = haoran_pan@qq.com
__date__ = 20210818
__version__ = v4.1
'''

import sys
import re
import pandas as pd


def countfamilies(file):
    df = pd.read_table(file,sep='\t')
    col_names = [col for col in df][1:-1]
    for colnames in col_names:
        species_name = re.findall('([0-9a-zA-Z\_\.\-]+)\.pep',colnames)[0]
        families_OG = df[df[colnames] !=0].iloc[:,0].tolist()
        with open(species_name+'.group.txt','w') as f:
            for OG in families_OG:
                f.write(OG+'\n')


if __name__ == '__main__':
    from optparse import OptionParser
    from rich.traceback import install
    install()
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args) == 1:
        count_file = args[0]
        countfamilies(count_file)
    else:
        sys.exit(p.print_help())
