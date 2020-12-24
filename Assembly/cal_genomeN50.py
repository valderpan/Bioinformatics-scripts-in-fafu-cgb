#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2020/12/24

import sys
sys.path.append('/public1/home/stu_panhaoran/scripts')
import fa2dict
import Fontcolor as Font

def calcultate_N50(genome_ctg_dict):
    ctgs_length = 0
    for key in genome_ctg_dict.keys():
        ctgs_length += len(genome_ctg_dict[key])
    Font.Log_output('ctg genome size:  {}'.format(ctgs_length))
    N50_length = round(ctgs_length/2)
    N50 = 0
    for key in sorted(genome_ctg_dict.items(),key=lambda x:len(x[1]),reverse=True):
        if N50 < N50_length:
            N50 += len(genome_ctg_dict[key[0]])
        else:
            Font.Log_output('N50 contig:  {} , length  {}'.format(key[0],len(genome_ctg_dict[key[0]])))
            break
if __name__ == '__main__':
    if len(sys.argv) != 2:
        Font.Scripts_tip('Calculate the N50 of the contig genome')
        print('Usage:')
        print('\tpython {0} <genome_ctg_file>'.format(sys.argv[0]))
    else:
        genomedict = fa2dict.Fa2dict(sys.argv[1])
        calcultate_N50(genomedict)
