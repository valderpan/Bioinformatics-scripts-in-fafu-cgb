#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/12/7


'''
%prog <RM out file> <genome_fasta>

Counting genomic repeat sequences predicted by RepeatMasker

>>> python %prog ZG.v20211208.genome.fasta.out ZG.v20211208.genome.fasta

__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20220402'
__version__ = 'v1.2'
'''


import sys
import fa2dict_Bio as fdB
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


def get_genomesize(fastaD):
    num = 0
    for key in fastaD.keys():
        num += len(fastaD[key])
    return num


def output_result(family_lenD,family_numD,genome_size):
    repeat_num = 0
    for key in family_lenD.keys():
        repeat_num += family_lenD[key]
    for key in family_lenD.keys():
        print('{}:\n\tNumber: {}\n\tLength (bp): {}\n\tPercent of genome(%): {}\n\tPecentage of repeat(%): {}'.format(
                                                           key,family_numD[key],family_lenD[key],
                                                           round(family_lenD[key]/int(genome_size)*100,3),
                                                           round(family_lenD[key]/repeat_num*100,3)))


if __name__ == '__main__':
    from optparse import OptionParser
    from rich.traceback import install
    install()
    p = OptionParser(__doc__)
    opts, args = p.parse_args()
    if len(args) == 2:
        out_file = args[0]
        genome_fasta = args[1]
        lenD, numD = parse_RMout(out_file)
        fastaD = fdB.Fa2dict_Bio(genome_fasta)
        genome_size = get_genomesize(fastaD)
        print('Genome size : {}bp'.format(genome_size))
        output_result(lenD,numD,genome_size)
    else:
        print('Usage:')
        print('\tpython {0} <RM.out>'.format(sys.argv[0]))
