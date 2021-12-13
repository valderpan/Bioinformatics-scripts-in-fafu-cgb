#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/12/7


'''
%prog <RM out file>

将RepeatMasker预测出来的基因组重复序列成立为表发级别的表格形式

__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20211213'
__version__ = 'v1.0'
'''


import sys
from rich import print


def parse_RMout(file):
    family_lenD = {};family_numD = {}
    with open(file,'r') as f:
        for line in f:
            if line.strip().startswith('SW') or line.strip().startswith('score'):
                continue
            else:
                line_list = line.strip().split()
                if line_list:
                    family = line_list[10]
                    start = line_list[5]
                    end = line_list[6]
                    length = int(end)-int(start)
                    if family not in family_lenD:
                        family_lenD[family] = 0
                        family_numD[family] = 0
                        family_lenD[family] += length
                        family_numD[family] += 1
                    else:
                        family_lenD[family] += length
                        family_numD[family] += 1
    return family_lenD,family_numD


def output_result(family_lenD,family_numD):
    for key in family_lenD.keys():
        print('{}:\n\tNumber: {}\n\tLength (bp): {}\n\tPercent (%): {}'.format(key,family_numD[key],
                                                           family_lenD[key],round(family_lenD[key]/6577007109*100,3)))


if __name__ == '__main__':
    from optparse import OptionParser
    from rich.traceback import install
    install()
    p = OptionParser(__doc__)
    opts, args = p.parse_args()
    if len(args) == 1:
        out_file = args[0]
        lenD, numD = parse_RMout(out_file)
        output_result(lenD, numD)
    else:
        print('Usage:')
        print('\tpython {0} <RM.out>'.format(sys.argv[0]))