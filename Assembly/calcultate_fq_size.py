#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2020/12/23

import Genome_phasing as phasing
import sys
import Fontcolor as Font
def calcultate_sample_size(file):
    phasing.calcultate_sample_length(file)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        Font.Scripts_tip('Calculate the amount of data for sequencing samples')
        print('Usage:')
        print('\tpython {0} <sample.fq.gz>'.format(sys.argv[0]))
    else:
        calcultate_sample_size(sys.argv[1])
