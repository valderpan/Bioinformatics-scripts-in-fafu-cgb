#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/3/3

#-----------------------------------------
#  CPCS: Contig position correct scripts |
#-----------------------------------------

import pandas as pd
import sys
import re

def remove_redundant_comment_line(group_file):
    col_names = ['ROC_id','ROC_chr','ROC_group','ROC_gene_position','Sb_id','Sb_chr','Sb_gene_position']
    chr_df = pd.read_excel(group_file, header=None,names=col_names)
    duplicate_row = []
    for index,row in chr_df.iterrows():
        if '#' in chr_df.iloc[index]['ROC_id'] and '#' in chr_df.iloc[index+1]['ROC_id']:
            duplicate_row.append(index+1)
    chr_df = chr_df.drop(labels=duplicate_row).reset_index(drop=True)
    return chr_df

def SeparateClassifiBlocks(group_df,query_list):
    split_row_list = []
    for index, row in group_df.iterrows():
        if group_df.loc[index]['ROC_id'] == '###':
            split_row_list.append(index)

    concat_row_list = []
    for i in range(0, len(split_row_list)):
        if i != len(split_row_list) - 1:
            group_freg = group_df.iloc[split_row_list[i]:split_row_list[i + 1], :]
            concat_row_list.append(group_freg)
        else:
            group_freg = group_df.iloc[split_row_list[i]:group_df.shape[0], :]
            concat_row_list.append(group_freg)

    ConGroup = str(list(query_list)[0]).split(',')
    print(ConGroup)
    nonconserverd_group_concat = []
    conserverd_group_concat = []
    for block in concat_row_list:
        group_id = re.findall('[0-9]+',block.iloc[1,:]['ROC_group'])[0]
        if group_id not in ConGroup:
            nonconserverd_group_concat.append(block)
        elif group_id in ConGroup:
            conserverd_group_concat.append(block)

    return ConGroup,nonconserverd_group_concat,conserverd_group_concat
