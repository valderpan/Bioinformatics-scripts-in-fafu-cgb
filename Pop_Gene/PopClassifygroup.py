#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/10/31

'''
%prog <pca.eigenvec> <pca num> <sample2group file> [pca.match.eigenvec]
This script is used to classify population samples

__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20201101'
__version__ = 'v1.1'
'''

import re
import sys
import pandas as pd


def Read_eigenvec(pca_file,pcas):
    pcanames = ['pc{}'.format(i) for i in range(1,int(pcas)+1)]
    pcadf = pd.read_table(pca_file,sep=' ',names=['SampleID','group']+pcanames)
    pca_prefix = re.findall('([0-9A-Za-z\_\-\.]+)\.eigenvec',pca_file)[0]
    return pcadf,pca_prefix


def Read_relation(file):
    df = pd.read_table(file,sep='\t',names=['SampleID','Group'])
    sample2group = df.set_index('SampleID').to_dict()['Group']
    return sample2group


def match_relation2eigenvec(pcadf,sample2group,pca_prefix):
    pcadf['group'] = pcadf['group'].map(sample2group)
    pcadf['group'] = pcadf['group'].map(str)
    pcadf.to_csv('{}.match.eigenvec'.format(pca_prefix),sep=' ',header=False,index=False)




if __name__ == '__main__':
    from optparse import OptionParser

    p = OptionParser(__doc__)
    opts, args = p.parse_args()
    if len(args) == 3:
        pca_file = args[0]
        pcas = args[1]
        relation_file = args[2]
        pcadf, pca_prefix = Read_eigenvec(pca_file,pcas)
        sample2group = Read_relation(relation_file)
        match_relation2eigenvec(pcadf, sample2group, pca_prefix)
    else:
        sys.exit(p.print_help())