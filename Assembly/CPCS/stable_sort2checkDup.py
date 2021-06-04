#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/6/1

'''
此脚本用来将stable中每个区块每个contig进行sort然后对stable内部进行去冗余,方便于解麻子
'''

import pandas as pd
import numpy as np
import re
import sys
sys.path.append(r'D:\pycharm\panhr\FAFU_CGB\Genomic_exercise')
import SetLog

def read_stable_file(stable_xlsx):
    df = pd.read_excel(stable_xlsx,names=['V1','V2','V3','V4','V5','V6','V7'])
    tmp_row = {'V1': "###", 'V2': np.nan, 'V3': np.nan, 'V4': np.nan, 'V5': np.nan, 'V6': np.nan, 'V7': np.nan}
    tmp_df = pd.DataFrame(tmp_row, pd.Index(range(1)))
    df = pd.concat([tmp_df, df], ignore_index=True)
    # df = df[~df['V1'].str.contains('###')]
    df_2stat = df.dropna(how='any')
    # print(df)
    df_name = re.findall('(C[0-9]+)\.stable*', stable_xlsx)[0]
    # print(df_name)
    return df,df_name,df_2stat.shape[0]


def stat_block_pos(stable_df):
    split_row_list = []
    for index, row in stable_df.iterrows():
        if stable_df.loc[index]['V1'] == '###':
            split_row_list.append(index)
    # print(len(split_row_list))
    return split_row_list

def split_every_block(stable_df,split_row_list):
    concat_row_list = []
    for i in range(0, len(split_row_list)):
        if i != len(split_row_list) - 1:
            group_freg = stable_df.iloc[split_row_list[i]:split_row_list[i + 1], :]
            concat_row_list.append(group_freg)
        else:
            group_freg = stable_df.iloc[split_row_list[i]:stable_df.shape[0], :]
            concat_row_list.append(group_freg)
    # print(len(concat_row_list))
    return concat_row_list

def sort_perBlock_perContig(concat_row_list,originStable_rows,output_prefix):
    result_concat = []
    for block in concat_row_list:
        contig2min = {}
        contigID = block['V2'].unique().tolist()
        for ID in contigID:
            contig2min[ID]=block[block['V2'] == ID]['V7'].min()
        tmp_row = {'V1': "###", 'V2': np.nan, 'V3': np.nan, 'V4': np.nan, 'V5': np.nan, 'V6': np.nan, 'V7': np.nan}
        tmp_df = pd.DataFrame(tmp_row, pd.Index(range(1)))
        result_concat.append(tmp_df)
        for key in sorted(contig2min.items(),key=lambda x:x[1],reverse=False):
            if key[0] != np.nan:
                query_contig_df = block[block['V2'] == key[0]]
                query_contig_sorted_df = query_contig_df.sort_values(by='V7')
                result_concat.append(query_contig_sorted_df)


    result = pd.concat(result_concat,ignore_index=True)
    result_2stat = result.dropna(how='any')
    print(result)
    if result_2stat.shape[0] == originStable_rows:
        print('Origin stable rows : {}\nSorted stable rows : {}'.format(originStable_rows, result_2stat.shape[0]))
        SetLog.info_out('Origin stable rows == Sorted stable rows | [check OK ~]')
        result.to_excel('{}.stable.sort.xlsx'.format(output_prefix),header=False,index=False)
    else:
        print('Origin stable rows : {}\nSorted stable rows : {}'.format(originStable_rows, result_2stat.shape[0]))
        SetLog.error_out('Origin stable rows != Sorted stable rows !!')



if __name__ == '__main__':
    from optparse import OptionParser
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args)==1:
        origin_stable_xlsx = args[0]
        stable_df, df_name, stable_df_rows = read_stable_file(origin_stable_xlsx)
        split_row_list = stat_block_pos(stable_df)
        concat_row_list = split_every_block(stable_df,split_row_list)
        sort_perBlock_perContig(concat_row_list, stable_df_rows, df_name)
    else:
        sys.exit(p.print_help())