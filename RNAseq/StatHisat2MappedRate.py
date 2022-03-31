#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/03/31


'''
通过指定--summary-file参数得到的结果，计算hisat2运行后的各指标的比对情况

>>> hisat2 --summary-file summary -x <ht2-idx> {-1 <m1> -2 <m2>} [-S <sam>]
>>> python %prog summary

__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20220331'
__version__ = 'v1.0'

'''


import re
import sys
from rich import print


def read_result(file):
    rowL = []
    with open(file) as f:
        lines = (line.strip() for line in f)
        for line in lines:
            rowL.append(line)
    return rowL


def re_findends(patten,rowL):
    for i in rowL:
        if i.endswith(patten):
            return int(re.findall('([0-9]+)',i)[0])
        else:
            continue


def parse_result(rowL):
    readsL = []
    single_reads = re_findends('reads; of these:',rowL)
    uniq_reads1 = re_findends('aligned concordantly exactly 1 time',rowL)
    uniq_reads2 = re_findends('aligned discordantly 1 time',rowL)
    uniq_reads3 = re_findends('aligned exactly 1 time',rowL)
    multi_reads1 = re_findends('aligned concordantly >1 times',rowL)
    multi_reads2 = re_findends('aligned >1 times',rowL)
    return [single_reads,uniq_reads1,uniq_reads2,uniq_reads3,multi_reads1,multi_reads2]


def cal_rate(reads_num_list):
    Total_reads = reads_num_list[0]*2
    uniq_reads = reads_num_list[1]*2+reads_num_list[2]*2+reads_num_list[3]
    uniq_rate = uniq_reads/Total_reads*100
    multi_reads = reads_num_list[4]*2+reads_num_list[5]
    multi_rate = multi_reads/Total_reads*100
    return [Total_reads,uniq_reads,uniq_rate,multi_reads,multi_rate]


def print_output(logfile,resL):
    print('Sample file : {}'.format(logfile))
    print('    Total reads num : {}'.format(resL[0]))
    print('    Uniq reads num : {}'.format(resL[1]))
    print('    Uniq reads rate : {}%'.format(round(resL[2],2)))
    print('    Multi reads num : {}'.format(resL[3]))
    print('    Multi reads rate : {}%'.format(round(resL[4],2)))


if __name__ == "__main__":
    from optparse import OptionParser
    from rich.traceback import install
    install()
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args)==1:
        log_file = args[0]
        rowL = read_result(log_file)
        readsL = parse_result(rowL)
        resL = cal_rate(readsL)
        print_output(log_file,resL)
    else:
        sys.exit(p.print_help())