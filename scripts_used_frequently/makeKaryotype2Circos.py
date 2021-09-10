#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/9/10

'''
%prog <ori_file> <species abbreviation> <color_name>

Used to make Karyotype files to make circos diagrams
Note: The newly generated file will replace the original file !!!
'''

import sys
import pandas as pd


def read_orifile(file):
    df = pd.read_table(file,names=['seqid','end'])
    return df


def mutate_col(df,abb):
    df['chr'] = 'chr'
    df['type'] = '-'
    seqid = ['{}{}'.format(abb,i) for i in range(1,df.shape[0]+1)]
    df['ID'] = seqid
    df['start'] = 0
    df['colormap'] = df['ID']
    new_df = df.loc[:,['chr','type','seqid','ID','start','end','colormap']]
    return new_df


def output_result(file,new_df):
    new_df.to_csv(file,sep='\t',header=False,index=False)

def color_info(new_df,color):
    L = []
    for index,i in enumerate(new_df['seqid'].tolist()):
        if index != 0 and index % 5 != 0:
            if index < len(new_df['seqid'].tolist())-1:
                L.append('{}={},'.format(i,color))
            else:
                L.append('{}={}'.format(i, color))
        elif index == 0:
            L.append('{}={},'.format(i,color))
        elif index != 0 and index % 5 == 0:
            if index < len(new_df['seqid'].tolist())-1:
                L.append('{}={},\\\n'.format(i,color))
            else:
                L.append('{}={}'.format(i, color))

    print(''.join(L))

if __name__ == '__main__':
    from optparse import OptionParser

    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args) == 3:
        ori_file = args[0]
        abb = args[1]
        color = args[2]
        df = read_orifile(ori_file)
        new_df = mutate_col(df,abb)
        output_result(ori_file,new_df)
        color_info(new_df,color)
    else:
        sys.exit(p.print_help())