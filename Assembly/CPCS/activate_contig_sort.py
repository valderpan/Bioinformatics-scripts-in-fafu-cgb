#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/4/22

'''
对stable去除冗余后并将冗余的移动到activate后，对activate的contig进行排序，方便于染色体矫正
'''

import pandas as pd
import numpy as np
import re
import sys
sys.path.append(r'D:\pycharm\panhr\FAFU_CGB\Genomic_exercise')
import SetLog


def read_moved_acitvate_file(moved_activate_xlsx):
    df = pd.read_excel(moved_activate_xlsx, header=None, sheet_name=0,comment='#')
    df.rename(columns={0:'V1',1:'V2',2:'V3',3:'V4',4:'V5',5:'V6',6:'V7'},inplace=True)
    df_2stat = df.dropna(how='any')
    print(df_2stat)
    print(df)
    df_name = re.findall('(C[0-9]+)\.activate*', moved_activate_xlsx)[0]
    return df,df_name,df_2stat.shape[0]


def store_activate_contig(activate_df):
    activate_contig = activate_df['V2'].unique().tolist()
    contig2min = {}
    for contig in activate_contig:
        contig_df_min = activate_df[activate_df['V2'] == contig]['V7'].min()
        contig2min[contig] = contig_df_min
    return contig2min


def sort_activate_contig(activate_df,congtig2min_Dict,output_prefix,origin_activate_rows):
    concat_list = []
    tmp_row = {'V1': "###", 'V2': np.nan, 'V3': np.nan, 'V4': np.nan, 'V5': np.nan, 'V6': np.nan, 'V7': np.nan}
    tmp_df = pd.DataFrame(tmp_row, pd.Index(range(1)))
    for item in sorted(congtig2min_Dict.items(),key=lambda x:x[1],reverse=False):
        concat_list.append(tmp_df)
        # print(item)
        tmp_query_df = activate_df[activate_df['V2'] == item[0]]
        tmp_query_df_sort = tmp_query_df.sort_values(by='V7')
        concat_list.append(tmp_query_df_sort)
        # break
    result = pd.concat(concat_list,ignore_index=True)
    print(result)
    tmp_result = result[~result['V1'].str.contains('###')]
    print(tmp_result)
    # tmp_result = result[result['V2'].str.contains('tig')]
    # print(result)
    if origin_activate_rows == tmp_result.shape[0]:
        SetLog.info_out('Origin activate rows == Sorted activate rows | [check OK ~]')
        result.to_excel('{}.activate.sort.xlsx'.format(output_prefix),header=False,index=False)
    else:
        print('Origin activate rows:{}\nSorted activate rows:{}'.format(origin_activate_rows,tmp_result.shape[0]))
        SetLog.error_out('Origin activate rows != Sorted activate rows !!')
        result.to_excel('{}.activate.wrong.sort.xlsx'.format(output_prefix),header=False,index=False)



if __name__ == '__main__':
    # activate_df,df_name,activate_rows = read_moved_acitvate_file(r'E:\Badila\synteny_analysis\Badila染色体第二次校正\C5.activate.xlsx')
    activate_df,df_name,activate_rows = read_moved_acitvate_file(sys.argv[1])
    contig2min_Dict = store_activate_contig(activate_df)
    sort_activate_contig(activate_df,contig2min_Dict,df_name,activate_rows)
