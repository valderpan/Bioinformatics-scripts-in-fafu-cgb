#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/6/3

'''
此脚本用来将stable自身去冗余后对stable中每个区块每个contig按照contig为单位进行sort并输出,方便于解麻子
'''

import numpy as np
import pandas as pd
import sys
import stable_sort2checkDup
sys.path.append(r'D:\pycharm\panhr\FAFU_CGB\Genomic_exercise')
import SetLog

def sort_perContig(concat_row_list,originStable_rows,output_prefix):
    result_concat = []
    for block in concat_row_list:
        tmp_row = {'V1': "###", 'V2': np.nan, 'V3': np.nan, 'V4': np.nan, 'V5': np.nan, 'V6': np.nan, 'V7': np.nan}
        tmp_df = pd.DataFrame(tmp_row, pd.Index(range(1)))
        contig2min = {}
        contigID = block['V2'].unique().tolist()
        for ID in contigID:
            contig2min[ID] = block[block['V2'] == ID]['V7'].min()


        for key in sorted(contig2min.items(), key=lambda x: x[1], reverse=False):
            if key[0] != np.nan:
                query_contig_df = block[block['V2'] == key[0]]
                query_contig_sorted_df = query_contig_df.sort_values(by='V7')
                result_concat.append(tmp_df)
                result_concat.append(query_contig_sorted_df)

    result = pd.concat(result_concat, ignore_index=True)
    print(result)
    result_2stat = result.dropna(how='any')
    print(result_2stat)
    if result_2stat.shape[0] == originStable_rows:
        print('Origin stable rows : {}\nSorted stable rows : {}'.format(originStable_rows, result_2stat.shape[0]))
        SetLog.info_out('Origin stable rows == Sorted stable rows | [check OK ~]')
        result.to_excel('{}.stable.sort.AftermoveDup.xlsx'.format(output_prefix), header=False, index=False)
    else:
        print('Origin stable rows : {}\nSorted stable rows : {}'.format(originStable_rows,result_2stat.shape[0]))
        SetLog.error_out('Origin stable rows != Sorted stable rows !!')


if __name__ == '__main__':
    from optparse import OptionParser

    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args) == 1:
        origin_stable_xlsx = args[0]
        stable_df, df_name, stable_df_rows = stable_sort2checkDup.read_stable_file(origin_stable_xlsx)
        split_row_list = stable_sort2checkDup.stat_block_pos(stable_df)
        concat_row_list = stable_sort2checkDup.split_every_block(stable_df,split_row_list)
        sort_perContig(concat_row_list,stable_df_rows,df_name)
