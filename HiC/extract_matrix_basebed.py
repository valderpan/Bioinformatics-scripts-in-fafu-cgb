#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/8/16

'''
%prog <bed> <matrix> <chrID>

Matrix extraction based on bed files
'''
import re
import sys
import pandas as pd


def read_bed(bed):
    df = pd.read_table(bed,sep='\t',names=['chr','start','end','index'],float_precision='round_trip')
    return df


def read_matrix(matrix):
    pre = re.findall('([0-9a-zA-Z\-\_\.]+)\.matrix',matrix)[0]
    df = pd.read_table(matrix,sep='\t',names=['bin1','bin2','score'],float_precision='round_trip')
    return df,pre


def get_Qregion(bed_df,chrID,pre):
    q_df = bed_df[bed_df['chr'] == chrID]
    i_l = q_df['index'].tolist()
    i_s = min(i_l)
    i_e = max(i_l)
    q_df.to_csv('{}.{}.bed'.format(pre,chrID),sep='\t',header=False,index=False)
    return i_s,i_e


def get_Qmatrix(matrix_df,i_s,i_e):
    q_df = matrix_df[(matrix_df['bin1'].map(int) >= int(i_s)) & (matrix_df['bin1'].map(int) <= int(i_e))]
    q_df = q_df[(q_df['bin2'].map(int) >= int(i_s)) & (q_df['bin2'].map(int) <= int(i_e))]
    return q_df


def output(q_df,chrID,pre):
    q_df['score'] = round(q_df['score'],6)
    q_df.to_csv('{}.{}.matrix'.format(pre,chrID),sep='\t',header=False,index=False)


if __name__ == '__main__':
    from rich.console import Console
    from rich.traceback import install
    from optparse import OptionParser


    console = Console()
    install()
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args) == 3:
        bed_file = args[0]
        matrix_file = args[1]
        chrID = args[2]
        bed = read_bed(bed_file)
        matrix,pre = read_matrix(matrix_file)
        i_s,i_e = get_Qregion(bed,chrID,pre)
        q_df = get_Qmatrix(matrix,i_s,i_e)
        output(q_df,chrID,pre)
    else:
        sys.exit(p.print_help())