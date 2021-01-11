#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/1/11

import sys
import Genome_phasing
sys.path.append(r'D:\pycharm\panhr\FAFU_CGB\Genomic_exercise')
import Fontcolor as Font

def calculate_genome_size(genome_fa):
    query_genome_size = Genome_phasing.cal_genome_size(genome_fa)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        Font.Scripts_tip('Calculate the size of query genome')
        print('Usage:')
        print('\tpython {0} <genome_fa>'.format(sys.argv[0]))
    else:
        calculate_genome_size(sys.argv[1])