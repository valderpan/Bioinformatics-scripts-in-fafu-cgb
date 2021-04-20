#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/02/28

import sys
sys.path.append('/share/home/stu_panhaoran/scripts')
import Fontcolor as Font
import fa2dict_Bio

def cal_ctgN50(CtgGenome_dict):
    ctgs_length = 0
    contigs_num = 0
    for key in CtgGenome_dict.keys():
        contigs_num += 1
        ctgs_length += len(CtgGenome_dict[key])

    Font.Log_output('Total number of ctgs:  {}'.format(contigs_num))
    Font.Log_output('Total length of ctgs:  {} bp'.format(ctgs_length))
    Font.Log_output('Average contig size:  {} bp'.format(round(ctgs_length / contigs_num, 2)))

    N50_length = ctgs_length / 2
    N50 = 0

    Font.Log_output('Longest ctg sequence length:  {} bp'.format(
        len(list(enumerate(sorted(CtgGenome_dict.items(), key=lambda x: len(x[1]))))[-1][1][1])))
    Font.Log_output('Shortest ctg sequence length:  {} bp'.format(
        len(list(enumerate(sorted(CtgGenome_dict.items(), key=lambda x: len(x[1]))))[0][1][1])))

    for num,key in enumerate(sorted(CtgGenome_dict.items(), key=lambda x: len(x[1]), reverse=True)):
        N50 += len(CtgGenome_dict[key[0]])
        if N50 >= N50_length:
            Font.Log_output('N50 contig:  {} ({} sequences),length  {} bp'.format(key[0],num+1 ,len(CtgGenome_dict[key[0]])))
            break


if len(sys.argv) != 2:
    Font.Scripts_tip('Calculate the N50 of the contig genome')
    print('Usage:')
    print('\tpython {0} <genome_ctg_file>'.format(sys.argv[0]))
else:
    genomedict = fa2dict_Bio.Fa2dict_Bio(sys.argv[1])
    Font.Log_output('Results for {}'.format(sys.argv[1]))
    cal_ctgN50(genomedict)
