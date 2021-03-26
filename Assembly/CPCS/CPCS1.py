#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/3/3

#-----------------------------------------
#  CPCS: Contig position correct scripts |
#-----------------------------------------
'''
删去每个group中片段之间多余的"#",按照"#"将每个片段分割，将保守、非保守的分别划分为一个dataframe
'''

import pandas as pd
import sys
import re

def remove_redundant_comment_line(group_file):
    '''
    删去每个group中片段之间多余的"#"
    输入:未校正的染色体文件，eg:g2.xlsx
    '''
    col_names = ['ROC_id','ROC_chr','ROC_group','ROC_gene_position','Sb_id','Sb_chr','Sb_gene_position']
    chr_df = pd.read_excel(group_file, header=None,names=col_names)
    duplicate_row = []
    for index,row in chr_df.iterrows():
        if '#' in chr_df.iloc[index]['ROC_id'] and '#' in chr_df.iloc[index+1]['ROC_id']:
            duplicate_row.append(index+1)
    chr_df = chr_df.drop(labels=duplicate_row).reset_index(drop=True)
    return chr_df

def SeparateClassifiBlocks(group_df,query_list=0):
    '''
    按照"#"将每个片段分割，将保守、非保守的分别划分为一个dataframe
    输入文件:
    group_file:函数remove_redundant_comment_line的结果
    query_list:保留下来的10个group的id
    return:10组保守groups的片段,非保守groups的片段
    '''
    split_row_list = []
    for index, row in group_df.iterrows():
        if group_df.loc[index]['ROC_id'] == '###':
            split_row_list.append(index)
    # print(split_row_list)
    concat_row_list = []
    for i in range(0, len(split_row_list)):
        if i != len(split_row_list) - 1:
            group_freg = group_df.iloc[split_row_list[i]:split_row_list[i + 1], :]
            concat_row_list.append(group_freg)
        else:
            group_freg = group_df.iloc[split_row_list[i]:group_df.shape[0], :]
            concat_row_list.append(group_freg)
    # print(concat_row_list)
    # print(len(concat_row_list)) #590
    ConGroup= [1,11,12,19,22,24,59,73]
    nonconserverd_group_concat = []
    conserverd_group_concat = []
    for block in concat_row_list:
        group_id = re.findall('[0-9]+',block.iloc[1,:]['ROC_group'])[0]
        if int(group_id) not in ConGroup:
            nonconserverd_group_concat.append(block)
        elif int(group_id) in ConGroup:
            conserverd_group_concat.append(block)
    # nonresult = pd.concat(nonconserverd_group_concat)
    # conresult = pd.concat(conserverd_group_concat)
    # nonresult.to_excel('nonresult.xlsx',header=False,index=False)
    # conresult.to_excel('conresult.xlsx',header=False,index=False)

    # print(nonresult)
    # print(conresult)
    return ConGroup,nonconserverd_group_concat,conserverd_group_concat






if __name__ == '__main__':
    chr_df = remove_redundant_comment_line(r'E:\Badila\synteny_analysis\C03\C03.sort.xlsx')
    SeparateClassifiBlocks(chr_df)
