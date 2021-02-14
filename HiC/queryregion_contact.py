#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/1/23


import pandas as pd
import sys
sys.path.append('/public1/home/stu_panhaoran/scripts')
import Fontcolor as Font

'''
Find the HiC contact intensity of the query region according to the matrix file after correction
'''

def get_queryregion1(matirx_file,start,end):
    df = pd.read_table(matirx_file,header=None,sep='\t',names=['firstbin','secondbin','frequency'],float_precision='round_trip')
    df_start = df[df['firstbin'] == int(start)].index[0]
    df_end = df[df['firstbin'] == int(end)].index[-1]
    df = df.iloc[df_start:df_end+1,]
    print(df)
    return df

def get_queryregion2(matirx,start,end,output):
    concat_list = []
    for i in matirx['firstbin'].unique().tolist():
        dd = matirx[matirx['firstbin'] == i]
        dd_start = dd[dd['secondbin'] == int(start)].index[0]
        dd_end = dd[dd['secondbin'] == int(end)].index[0]
        df = dd.loc[dd_start:dd_end]
        concat_list.append(df)
    result = pd.concat(concat_list)
    print(result)
    result.to_csv(output,sep='\t',header=False,index=False)


if __name__ == '__main__':
    if len(sys.argv) != 7:
        Font.Scripts_tip('Get query region contact according to Hic-Pro contact matrix')
        print('Usage:')
        print('\tpython {0} <corrected.matrix> firstbin_start firstbin_end secondbin_start secondbin_end <output.matrix>'.format(sys.argv[0]))
    else:
        df = get_queryregion1(sys.argv[1],sys.argv[2],sys.argv[3])
        get_queryregion2(df,sys.argv[4],sys.argv[5],sys.argv[6])
