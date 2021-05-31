#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/5/30

'''
此脚本用于根据共线性图将[高粱每条染色体:甘蔗每条group]之间的anchors文件(eg:C1.anchors.sorted.xlsx)进行处理
输出文件:每个ref_chr的stable_xlsx与activate_xlsx,对这两个文件进行手动矫正
'''

import pandas as pd
import re
import numpy as np
import sys
sys.path.append(r'D:\pycharm\panhr\FAFU_CGB\Genomic_exercise')
import SetLog
import argparse



__author__ = 'Haoran Pan'
__mail__ = 'haoran_pan@qq.com'
__date__ = '20210531'
__version__ = 'v1.0'


def read_origin_file(anchors_sort_xlsx):
    df = pd.read_excel(anchors_sort_xlsx,names=['V1','V2','V3','V4','V5','V6','V7'])
    tmp_row = {'V1':"###",'V2':np.nan,'V3':np.nan,'V4':np.nan,'V5':np.nan,'V6':np.nan,'V7':np.nan}
    tmp_df = pd.DataFrame(tmp_row,pd.Index(range(1)))
    df = pd.concat([tmp_df,df],ignore_index=True)
    df_name = re.findall('(C[0-9]+)\.anchors*',anchors_sort_xlsx)[0]
    # print(df_name)
    return df,df_name


def get_stablegroupID(query_list):
    StableID = str(list(query_list)[0]).split(',')
    print(StableID)
    return  StableID


def StoreandClassfiBlocks(origin_df,stableID):
    '''
    按照"#"将每个片段分割，将保守、非保守的分别划分为一个dataframe
    输入文件:
    group_file:函数remove_redundant_comment_line的结果
    query_list:保留下来的10个group的id
    return:10组保守groups的片段,非保守groups的片段
    '''
    split_row_list = []
    for index, row in origin_df.iterrows():
        if origin_df.loc[index]['V1'] == '###':
            split_row_list.append(index)
    # print(split_row_list)
    concat_row_list = []
    for i in range(0, len(split_row_list)):
        if i != len(split_row_list) - 1:
            group_freg = origin_df.iloc[split_row_list[i]:split_row_list[i + 1], :]
            concat_row_list.append(group_freg)
        else:
            group_freg = origin_df.iloc[split_row_list[i]:origin_df.shape[0], :]
            concat_row_list.append(group_freg)
    # print(concat_row_list[0])
    # print(len(concat_row_list))
    activate_group_concat = []
    stable_group_concat = []
    for block in concat_row_list:
        group_id = re.findall('[0-9]+', block.iloc[1, :]['V3'])[0]
        if group_id not in stableID:
            activate_group_concat.append(block)
        elif group_id in stableID:
            stable_group_concat.append(block)
    # print(activate_group_concat)
    # print(len(activate_group_concat))
    # print(len(stable_group_concat))
    return activate_group_concat,stable_group_concat



def Output_stable(stable_group_list):
    stable_df = pd.concat(stable_group_list)
    # print(stable_df)
    return stable_df


def OutputandSort_activate(activate_group_list):
    block2min = {}
    for i,j in enumerate(activate_group_list):
        min_value = j['V7'].min()
        # print(min_value)
        block2min[i] = min_value

    concat_list =[]
    for key in sorted(block2min.items(),key=lambda x:x[1],reverse=False):
        # print(key)
        # print(block2min[key[0]])
        # break
        concat_list.append(activate_group_list[key[0]])
    # print(len(concat_list))
    if len(activate_group_list) != len(concat_list):
        SetLog.error_out('Activate block number != Activate block number after storing in dictionary')
    elif len(activate_group_list) == len(concat_list):
        result = pd.concat(concat_list)
        # print(result)
        return result


def Output_result(df_name,stable_result,activate_result,output_path=0):
    if output_path:
        stable_result.to_excel('{}/{}.stable.xlsx'.format(output_path,df_name), header=False, index=False)
        activate_result.to_excel('{}/{}.activate.xlsx'.format(output_path,df_name), header=False, index=False)
    else:
        stable_result.to_excel('{}.stable.xlsx'.format(df_name),header=False,index=False)
        activate_result.to_excel('{}.activate.xlsx'.format(df_name),header=False,index=False)
    # pass


def main(args):
    input_origin_file = args.input
    stableID_list = args.query
    output_path = args.output
    df, df_name = read_origin_file(input_origin_file)
    stableID = get_stablegroupID(stableID_list)
    activate_group_list, stable_group_list = StoreandClassfiBlocks(df, stableID)
    stable_result = Output_stable(stable_group_list)
    activate_result = OutputandSort_activate(activate_group_list)
    Output_result(df_name, stable_result, activate_result, output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''This script is used to process the anchors file (eg:C1.anchors.sorted.xlsx) between [sorghum per chromosome:sugarcane per group] based on the synteny dotplot''',
        usage="python {} -i <anchors.sorted.tab> -q stableID -o [output_path] ".format(sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__, __mail__, __date__,
                                                                            __version__))
    parser.add_argument('-i', '--input', required=True, help='Input anchors.sort.tab')
    parser.add_argument('-q', '--query',required=True,nargs='+' ,help='Input stable group ID')
    parser.add_argument('-o', '--output',help='Specify output file path')
    args = parser.parse_args()
    main(args)


# if __name__ == '__main__':
#     df,df_name = read_origin_file(r'E:\Badila\synteny_analysis\Badila染色体第二次校正\C2.anchors.sorted.xlsx')
#     stableID = get_stablegroupID('9','10','11','12','13','14','15','16')
#     activate_group_list,stable_group_list = StoreandClassfiBlocks(df,stableID)
#     stable_result = Output_stable(stable_group_list)
#     activate_result = OutputandSort_activate(activate_group_list)
#     Output_result(df_name,stable_result,activate_result,r'E:\Badila\synteny_analysis\Badila染色体第二次校正')