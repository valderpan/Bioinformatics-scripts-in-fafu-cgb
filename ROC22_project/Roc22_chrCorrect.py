#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# @Time : 2020/9/18 10:02


import pandas as pd
import numpy as np
from intervaltree import Interval, IntervalTree
import argparse
import re
import xlrd


__author__ = 'Haoran Pan'
__mail__ = 'haoran_pan@qq.com'
__version__ = 'v1.0'

def remove_redundant_comment_line(group_file):
    chr_df = pd.read_excel(group_file,header=None)
    chr_df.rename(columns = {0:'ROC_id',1:'ROC_chr',2:'ROC_group',3:'ROC_gene_position',
                         4:'Sb_id',5:'Sb_chr',6:'Sb_gene_position'},inplace=True)
    redundant_row_list = []
    for index,row in chr_df.iterrows():
        if chr_df.loc[index]['ROC_id'] == '###':
            if chr_df.loc[index+1]['ROC_id'] == '###':
                redundant_row_list.append(index+1)
    chr_df.drop(labels=redundant_row_list,inplace=True)
    return chr_df

def Separate_groupAndSort_blocks(group_file,query_list):
    group_df = remove_redundant_comment_line(group_file)
    group_df = group_df.reset_index(drop=True)
    # print(group_df)
    split_row_list = []
    for index,row in group_df.iterrows():
        if group_df.loc[index]['ROC_id'] == '###':
            split_row_list.append(index)

    concat_row_list = []
    for i in range(0,len(split_row_list)):
        if i != len(split_row_list)-1:
            group_freg = group_df.iloc[split_row_list[i]:split_row_list[i+1],:]
            # print(group_freg)
            # group_freg.sort_values(by='Sb_gene_position',inplace=True) #同一片段里基因位置不能改动
            concat_row_list.append(group_freg)
        else:
            group_freg = group_df.iloc[split_row_list[i]:group_df.shape[0],:]
            # print(group_freg)
            # group_freg.sort_values(by='Sb_gene_position', inplace=True) #同一片段里基因位置不能改动
            concat_row_list.append(group_freg)

    tmp_query_list = list(query_list)
    conserverd_group_num =tmp_query_list[0].split(',')
    # print(conserverd_group_num)
    nonconserverd_group_concat = []
    conserverd_group_concat = []
    for i in concat_row_list:
        group_id = re.findall('[0-9]+',i.iloc[1,:]['ROC_group'])[0]
        if group_id not in conserverd_group_num:
            nonconserverd_group_concat.append(i)
        else:
            conserverd_group_concat.append(i)
    nonconserverd_group = pd.concat(nonconserverd_group_concat)
    conserved_group = pd.concat(conserverd_group_concat)

    # --------------------将10组保守的group排序后输出--------------------------------------
    writer_conserverd = pd.ExcelWriter('conserved_group.xlsx')
    for num in conserverd_group_num:
        # print(num)
        conserverd_group_individual_concat = []
        for j in conserverd_group_concat:
            group_id = re.findall('[0-9]+',j.iloc[1,:]['ROC_group'])[0]
            if group_id == num:
                # print(j)
                conserverd_group_individual_concat.append(j)
        conserverd_group_individual=pd.concat(conserverd_group_individual_concat)
        # print(conserverd_group_individual)  #将保守的10个group的dataframe进行输出
        sub_split_blocks_row = []
        for index, row in conserverd_group_individual.iterrows():
            if group_df.loc[index]['ROC_id'] == '###':
                sub_split_blocks_row.append(index)

        group_split_blocks_dict = {}
        for i in range(0, len(sub_split_blocks_row)):
            if i != len(sub_split_blocks_row) -1:
                group_split_blocks_dict[i] = group_df.iloc[sub_split_blocks_row[i]:sub_split_blocks_row[i + 1], :]
            else:
                group_split_blocks_dict[i] = group_df.iloc[sub_split_blocks_row[i]:conserverd_group_individual.index[-1]+1, :]

        tmp = {}
        for j in group_split_blocks_dict.keys():
            tmp[int(group_split_blocks_dict[j]['Sb_gene_position'].mean())] = group_split_blocks_dict[j]

        sorted_blocks_concat = []
        for i in sorted(tmp):
            sorted_blocks_concat.append(tmp[i])
        sorted_blocks = pd.concat(sorted_blocks_concat)

        sorted_blocks.to_excel(writer_conserverd,sheet_name=str(num),header=None,index=False)
    writer_conserverd.save()
    #---------------------------------------------------------------------------------

    #-------------------------------输出非保守的每个group，每个excel包含10个group，每个group占一个sheet------------------------
    ##------1-20-------
    writer_nonconserverd_1_20 = pd.ExcelWriter('nonconserved_group_1-20.xlsx')
    for num in range(1,21):
        if str(num) not in conserverd_group_num:
            nonconserverd_group_individual_concat = []
            for j in nonconserverd_group_concat:
                group_id = re.findall('[0-9]+', j.iloc[1, :]['ROC_group'])[0]
                if group_id == str(num):
                    nonconserverd_group_individual_concat.append(j)
            try:
                nonconserverd_group_individual = pd.concat(nonconserverd_group_individual_concat)
            except:
                continue
            sub_split_blocks_row = []
            for index, row in nonconserverd_group_individual.iterrows():
                if group_df.loc[index]['ROC_id'] == '###':
                    sub_split_blocks_row.append(index)

            group_split_blocks_dict = {}
            for i in range(0, len(sub_split_blocks_row)):
                if i != len(sub_split_blocks_row) -1:
                    group_split_blocks_dict[i] = group_df.iloc[sub_split_blocks_row[i]:sub_split_blocks_row[i + 1], :]
                else:
                    group_split_blocks_dict[i] = group_df.iloc[sub_split_blocks_row[i]:nonconserverd_group_individual.index[-1]+1, :]
            tmp = {}
            for j in group_split_blocks_dict.keys():
                tmp[int(group_split_blocks_dict[j]['Sb_gene_position'].mean())] = group_split_blocks_dict[j]

            sorted_blocks_concat = []
            for i in sorted(tmp):
                sorted_blocks_concat.append(tmp[i])
            sorted_blocks = pd.concat(sorted_blocks_concat)
            sorted_blocks.to_excel(writer_nonconserverd_1_20, sheet_name=str(num), header=None, index=False)
    writer_nonconserverd_1_20.save()
    #------21-40-------
    writer_nonconserverd_21_40 = pd.ExcelWriter('nonconserved_group_21-40.xlsx')
    for num in range(21, 41):
        if str(num) not in conserverd_group_num:
            nonconserverd_group_individual_concat = []
            for j in nonconserverd_group_concat:
                group_id = re.findall('[0-9]+', j.iloc[1, :]['ROC_group'])[0]
                if group_id == str(num):
                    nonconserverd_group_individual_concat.append(j)
            try:
                nonconserverd_group_individual = pd.concat(nonconserverd_group_individual_concat)
            except:
                continue
            sub_split_blocks_row = []
            for index, row in nonconserverd_group_individual.iterrows():
                if group_df.loc[index]['ROC_id'] == '###':
                    sub_split_blocks_row.append(index)

            group_split_blocks_dict = {}
            for i in range(0, len(sub_split_blocks_row)):
                if i != len(sub_split_blocks_row) - 1:
                    group_split_blocks_dict[i] = group_df.iloc[sub_split_blocks_row[i]:sub_split_blocks_row[i + 1], :]
                else:
                    group_split_blocks_dict[i] = group_df.iloc[sub_split_blocks_row[i]:nonconserverd_group_individual.index[-1] + 1,:]
            tmp = {}
            for j in group_split_blocks_dict.keys():
                tmp[int(group_split_blocks_dict[j]['Sb_gene_position'].mean())] = group_split_blocks_dict[j]

            sorted_blocks_concat = []
            for i in sorted(tmp):
                sorted_blocks_concat.append(tmp[i])
            sorted_blocks = pd.concat(sorted_blocks_concat)
            sorted_blocks.to_excel(writer_nonconserverd_21_40, sheet_name=str(num), header=None, index=False)
    writer_nonconserverd_21_40.save()
    # ------41-60-------
    writer_nonconserverd_41_60 = pd.ExcelWriter('nonconserved_group_41-60.xlsx')
    for num in range(41, 61):
        if str(num) not in conserverd_group_num:
            nonconserverd_group_individual_concat = []
            for j in nonconserverd_group_concat:
                group_id = re.findall('[0-9]+', j.iloc[1, :]['ROC_group'])[0]
                if group_id == str(num):
                    nonconserverd_group_individual_concat.append(j)
            try:
                nonconserverd_group_individual = pd.concat(nonconserverd_group_individual_concat)
            except:
                continue
            sub_split_blocks_row = []
            for index, row in nonconserverd_group_individual.iterrows():
                if group_df.loc[index]['ROC_id'] == '###':
                    sub_split_blocks_row.append(index)

            group_split_blocks_dict = {}
            for i in range(0, len(sub_split_blocks_row)):
                if i != len(sub_split_blocks_row) - 1:
                    group_split_blocks_dict[i] = group_df.iloc[sub_split_blocks_row[i]:sub_split_blocks_row[i + 1], :]
                else:
                    group_split_blocks_dict[i] = group_df.iloc[sub_split_blocks_row[i]:nonconserverd_group_individual.index[-1] + 1,:]
            tmp = {}
            for j in group_split_blocks_dict.keys():
                tmp[int(group_split_blocks_dict[j]['Sb_gene_position'].mean())] = group_split_blocks_dict[j]

            sorted_blocks_concat = []
            for i in sorted(tmp):
                sorted_blocks_concat.append(tmp[i])
            sorted_blocks = pd.concat(sorted_blocks_concat)
            sorted_blocks.to_excel(writer_nonconserverd_41_60, sheet_name=str(num), header=None, index=False)
    writer_nonconserverd_41_60.save()
    # ------61-80-------
    writer_nonconserverd_61_80 = pd.ExcelWriter('nonconserved_group_61-80.xlsx')
    for num in range(61, 81):
        if str(num) not in conserverd_group_num:
            nonconserverd_group_individual_concat = []
            for j in nonconserverd_group_concat:
                group_id = re.findall('[0-9]+', j.iloc[1, :]['ROC_group'])[0]
                if group_id == str(num):
                    nonconserverd_group_individual_concat.append(j)
            try:
                nonconserverd_group_individual = pd.concat(nonconserverd_group_individual_concat)
            except:
                continue
            sub_split_blocks_row = []
            for index, row in nonconserverd_group_individual.iterrows():
                if group_df.loc[index]['ROC_id'] == '###':
                    sub_split_blocks_row.append(index)

            group_split_blocks_dict = {}
            for i in range(0, len(sub_split_blocks_row)):
                if i != len(sub_split_blocks_row) - 1:
                    group_split_blocks_dict[i] = group_df.iloc[sub_split_blocks_row[i]:sub_split_blocks_row[i + 1], :]
                else:
                    group_split_blocks_dict[i] = group_df.iloc[sub_split_blocks_row[i]:nonconserverd_group_individual.index[-1] + 1,:]
            tmp = {}
            for j in group_split_blocks_dict.keys():
                tmp[int(group_split_blocks_dict[j]['Sb_gene_position'].mean())] = group_split_blocks_dict[j]

            sorted_blocks_concat = []
            for i in sorted(tmp):
                sorted_blocks_concat.append(tmp[i])
            sorted_blocks = pd.concat(sorted_blocks_concat)
            sorted_blocks.to_excel(writer_nonconserverd_61_80, sheet_name=str(num), header=None, index=False)
    writer_nonconserverd_61_80.save()
    # ------81-107-------
    writer_nonconserverd_81_end = pd.ExcelWriter('nonconserved_group_81-end.xlsx')
    for num in range(81, 110):
        if str(num) not in conserverd_group_num:
            nonconserverd_group_individual_concat = []
            for j in nonconserverd_group_concat:
                group_id = re.findall('[0-9]+', j.iloc[1, :]['ROC_group'])[0]
                if group_id == str(num):
                    nonconserverd_group_individual_concat.append(j)
            try:
                nonconserverd_group_individual = pd.concat(nonconserverd_group_individual_concat)
            except:
                continue
            sub_split_blocks_row = []
            for index, row in nonconserverd_group_individual.iterrows():
                if group_df.loc[index]['ROC_id'] == '###':
                    sub_split_blocks_row.append(index)

            group_split_blocks_dict = {}
            for i in range(0, len(sub_split_blocks_row)):
                if i != len(sub_split_blocks_row) - 1:
                    group_split_blocks_dict[i] = group_df.iloc[sub_split_blocks_row[i]:sub_split_blocks_row[i + 1], :]
                else:
                    group_split_blocks_dict[i] = group_df.iloc[
                                                 sub_split_blocks_row[i]:nonconserverd_group_individual.index[-1] + 1,
                                                 :]
            tmp = {}
            for j in group_split_blocks_dict.keys():
                tmp[int(group_split_blocks_dict[j]['Sb_gene_position'].mean())] = group_split_blocks_dict[j]

            sorted_blocks_concat = []
            for i in sorted(tmp):
                sorted_blocks_concat.append(tmp[i])
            sorted_blocks = pd.concat(sorted_blocks_concat)
            sorted_blocks.to_excel(writer_nonconserverd_81_end, sheet_name=str(num), header=None, index=False)
    writer_nonconserverd_81_end.save()


def correct_chr(conserverd_group,nonconserverd_group,nonconserverd_sheet,output_file):

    #-------------------------------------------------------
    ###-----nonconserverd----------------
    nonconserverd_df = pd.read_excel(nonconserverd_group, sheet_name=int(nonconserverd_sheet), names=['ROC_id','ROC_chr','ROC_group','ROC_gene_position','Sb_id',
                                                                               'Sb_chr','Sb_gene_position'],header=None)
    nonconserverd_split_blocks_row = []
    for index,row in nonconserverd_df.iterrows():
        if nonconserverd_df.loc[index]['ROC_id'] == '###':
            nonconserverd_split_blocks_row.append(index)
    non_group_id = re.findall('[0-9]+', nonconserverd_df.iloc[1, :]['ROC_group'])[0]
    nonconserverd_split_blocks_dict = {}
    for i in range(0,len(nonconserverd_split_blocks_row)):
        if i != len(nonconserverd_split_blocks_row)-1:
            nonconserverd_split_blocks_dict[i] = nonconserverd_df.iloc[nonconserverd_split_blocks_row[i]:nonconserverd_split_blocks_row[i+1],:]
        else:
            nonconserverd_split_blocks_dict[i] = nonconserverd_df.iloc[nonconserverd_split_blocks_row[i]:nonconserverd_df.shape[0],:]
    # print(nonconserverd_split_blocks_dict)

    # -------------将非保守的group按照max、min值构建区间-------------
    nonconserverd_blocks_interval_dict = {}
    for i in nonconserverd_split_blocks_dict.keys():
        block_start = nonconserverd_split_blocks_dict[i]['Sb_gene_position'].min()
        block_end   = nonconserverd_split_blocks_dict[i]['Sb_gene_position'].max()
        nonconserverd_blocks_interval_dict[Interval(block_start,block_end)] = nonconserverd_split_blocks_dict[i]
    # print(nonconserverd_blocks_interval_dict)
    #----------------------------------

    writer = pd.ExcelWriter(output_file)
    #---------conserverd1----------------
    conserverd_df_1 = pd.read_excel(conserverd_group, sheet_name=0, names=['ROC_id','ROC_chr','ROC_group','ROC_gene_position','Sb_id',
                                                                               'Sb_chr','Sb_gene_position'],header=None)
    #----------拆分为每一块---------------
    conserverd_split_blocks_row_1 = []
    for index,row in conserverd_df_1.iterrows():
        if conserverd_df_1.loc[index]['ROC_id'] == '###':
            conserverd_split_blocks_row_1.append(index)
    conserverd_split_blocks_dict_1 = {}
    for i in range(0,len(conserverd_split_blocks_row_1)):
        if i != len(conserverd_split_blocks_row_1)-1:
            conserverd_split_blocks_dict_1[i] = conserverd_df_1.iloc[conserverd_split_blocks_row_1[i]:conserverd_split_blocks_row_1[i+1],:]
        else:
            conserverd_split_blocks_dict_1[i] = conserverd_df_1.iloc[conserverd_split_blocks_row_1[i]:conserverd_df_1.shape[0],:]
    # print(conserverd_split_blocks_dict)
    #-------------将保守的group按照max、min值构建区间-------------
    conserverd_blocks_interval_dict_1 = {}
    for i in conserverd_split_blocks_dict_1.keys():
        block_start = conserverd_split_blocks_dict_1[i]['Sb_gene_position'].min()
        block_end   = conserverd_split_blocks_dict_1[i]['Sb_gene_position'].max()
        conserverd_blocks_interval_dict_1[Interval(block_start,block_end)] = conserverd_split_blocks_dict_1[i]
    # print(conserverd_blocks_interval_dict_1)

    #----------------区间判断----------------------------
    conserverd_blocks_interval_intertree_1 = IntervalTree(conserverd_blocks_interval_dict_1)
    not_overlaps_inter_1 = []
    for inter in nonconserverd_blocks_interval_dict.keys():
        if not conserverd_blocks_interval_intertree_1.overlaps(inter):
            not_overlaps_inter_1.append(inter)
    # print(not_overlaps_inter_1)
    ##----把需要转移的block进行转移------------------------
    new_moved_dict_1 = {}
    for i in not_overlaps_inter_1:
        new_moved_dict_1[i] = nonconserverd_blocks_interval_dict[i]
        conserverd_blocks_interval_dict_1.update(new_moved_dict_1)
        nonconserverd_blocks_interval_dict.pop(i)
    tmp1 = []
    for i in sorted(conserverd_blocks_interval_dict_1.keys()):
        tmp1.append(conserverd_blocks_interval_dict_1[i])
    result1 = pd.concat(tmp1)
    # print(result1)
    result1.to_excel(writer,sheet_name='1',header=None,index=False)

    #---------conserverd2----------------
    conserverd_df_2 = pd.read_excel(conserverd_group, sheet_name=1, names=['ROC_id','ROC_chr','ROC_group','ROC_gene_position','Sb_id',
                                                                               'Sb_chr','Sb_gene_position'],header=None)
    #----------拆分为每一块---------------
    conserverd_split_blocks_row_2 = []
    for index,row in conserverd_df_2.iterrows():
        if conserverd_df_2.loc[index]['ROC_id'] == '###':
            conserverd_split_blocks_row_2.append(index)
    conserverd_split_blocks_dict_2 = {}
    for i in range(0,len(conserverd_split_blocks_row_2)):
        if i != len(conserverd_split_blocks_row_2)-1:
            conserverd_split_blocks_dict_2[i] = conserverd_df_2.iloc[conserverd_split_blocks_row_2[i]:conserverd_split_blocks_row_2[i+1],:]
        else:
            conserverd_split_blocks_dict_2[i] = conserverd_df_2.iloc[conserverd_split_blocks_row_2[i]:conserverd_df_2.shape[0],:]
    # print(conserverd_split_blocks_dict)
    #-------------将保守的group按照max、min值构建区间-------------
    conserverd_blocks_interval_dict_2 = {}
    for i in conserverd_split_blocks_dict_2.keys():
        block_start = conserverd_split_blocks_dict_2[i]['Sb_gene_position'].min()
        block_end   = conserverd_split_blocks_dict_2[i]['Sb_gene_position'].max()
        conserverd_blocks_interval_dict_2[Interval(block_start,block_end)] = conserverd_split_blocks_dict_2[i]
    # print(conserverd_blocks_interval_dict_2)

    #----------------区间判断----------------------------
    conserverd_blocks_interval_intertree_2 = IntervalTree(conserverd_blocks_interval_dict_2)
    not_overlaps_inter_2 = []
    for inter in nonconserverd_blocks_interval_dict.keys():
        if not conserverd_blocks_interval_intertree_2.overlaps(inter):
            not_overlaps_inter_2.append(inter)
    # print(not_overlaps_inter_2)
    ##----把需要转移的block进行转移------------------------
    new_moved_dict_2 = {}
    for i in not_overlaps_inter_2:
        new_moved_dict_2[i] = nonconserverd_blocks_interval_dict[i]
        conserverd_blocks_interval_dict_2.update(new_moved_dict_2)
        nonconserverd_blocks_interval_dict.pop(i)
    tmp2 = []
    for i in sorted(conserverd_blocks_interval_dict_2.keys()):
        tmp2.append(conserverd_blocks_interval_dict_2[i])
    result2 = pd.concat(tmp2)
    # prin(result2)
    result2.to_excel(writer, sheet_name='2', header=None, index=False)


    #---------conserverd3----------------
    conserverd_df_3 = pd.read_excel(conserverd_group, sheet_name=2, names=['ROC_id','ROC_chr','ROC_group','ROC_gene_position','Sb_id',
                                                                               'Sb_chr','Sb_gene_position'],header=None)
    #----------拆分为每一块---------------
    conserverd_split_blocks_row_3 = []
    for index,row in conserverd_df_3.iterrows():
        if conserverd_df_3.loc[index]['ROC_id'] == '###':
            conserverd_split_blocks_row_3.append(index)
    conserverd_split_blocks_dict_3 = {}
    for i in range(0,len(conserverd_split_blocks_row_3)):
        if i != len(conserverd_split_blocks_row_3)-1:
            conserverd_split_blocks_dict_3[i] = conserverd_df_3.iloc[conserverd_split_blocks_row_3[i]:conserverd_split_blocks_row_3[i+1],:]
        else:
            conserverd_split_blocks_dict_3[i] = conserverd_df_3.iloc[conserverd_split_blocks_row_3[i]:conserverd_df_3.shape[0],:]
    # print(conserverd_split_blocks_dict)
    #-------------将保守的group按照max、min值构建区间-------------
    conserverd_blocks_interval_dict_3 = {}
    for i in conserverd_split_blocks_dict_3.keys():
        block_start = conserverd_split_blocks_dict_3[i]['Sb_gene_position'].min()
        block_end   = conserverd_split_blocks_dict_3[i]['Sb_gene_position'].max()
        conserverd_blocks_interval_dict_3[Interval(block_start,block_end)] = conserverd_split_blocks_dict_3[i]
    # print(conserverd_blocks_interval_dict_3)

    #----------------区间判断----------------------------
    conserverd_blocks_interval_intertree_3 = IntervalTree(conserverd_blocks_interval_dict_3)
    not_overlaps_inter_3 = []
    for inter in nonconserverd_blocks_interval_dict.keys():
        if not conserverd_blocks_interval_intertree_3.overlaps(inter):
            not_overlaps_inter_3.append(inter)
    # print(not_overlaps_inter_3)
    ##----把需要转移的block进行转移------------------------
    new_moved_dict_3 = {}
    for i in not_overlaps_inter_3:
        new_moved_dict_3[i] = nonconserverd_blocks_interval_dict[i]
        conserverd_blocks_interval_dict_3.update(new_moved_dict_3)
        nonconserverd_blocks_interval_dict.pop(i)
    tmp3 = []
    for i in sorted(conserverd_blocks_interval_dict_3.keys()):
        tmp3.append(conserverd_blocks_interval_dict_3[i])
    result3 = pd.concat(tmp3)
    # print(result3)
    result3.to_excel(writer, sheet_name='3', header=None, index=False)



    #---------conserverd4----------------
    conserverd_df_4 = pd.read_excel(conserverd_group, sheet_name=3, names=['ROC_id','ROC_chr','ROC_group','ROC_gene_position','Sb_id',
                                                                               'Sb_chr','Sb_gene_position'],header=None)
    #----------拆分为每一块---------------
    conserverd_split_blocks_row_4 = []
    for index,row in conserverd_df_4.iterrows():
        if conserverd_df_4.loc[index]['ROC_id'] == '###':
            conserverd_split_blocks_row_4.append(index)
    conserverd_split_blocks_dict_4 = {}
    for i in range(0,len(conserverd_split_blocks_row_4)):
        if i != len(conserverd_split_blocks_row_4)-1:
            conserverd_split_blocks_dict_4[i] = conserverd_df_4.iloc[conserverd_split_blocks_row_4[i]:conserverd_split_blocks_row_4[i+1],:]
        else:
            conserverd_split_blocks_dict_4[i] = conserverd_df_4.iloc[conserverd_split_blocks_row_4[i]:conserverd_df_4.shape[0],:]
    # print(conserverd_split_blocks_dict)
    #-------------将保守的group按照max、min值构建区间-------------
    conserverd_blocks_interval_dict_4 = {}
    for i in conserverd_split_blocks_dict_4.keys():
        block_start = conserverd_split_blocks_dict_4[i]['Sb_gene_position'].min()
        block_end   = conserverd_split_blocks_dict_4[i]['Sb_gene_position'].max()
        conserverd_blocks_interval_dict_4[Interval(block_start,block_end)] = conserverd_split_blocks_dict_4[i]
    # print(conserverd_blocks_interval_dict_4)

    #----------------区间判断----------------------------
    conserverd_blocks_interval_intertree_4 = IntervalTree(conserverd_blocks_interval_dict_4)
    not_overlaps_inter_4 = []
    for inter in nonconserverd_blocks_interval_dict.keys():
        if not conserverd_blocks_interval_intertree_4.overlaps(inter):
            not_overlaps_inter_4.append(inter)
    # print(not_overlaps_inter_4)
    ##----把需要转移的block进行转移------------------------
    new_moved_dict_4 = {}
    for i in not_overlaps_inter_4:
        new_moved_dict_4[i] = nonconserverd_blocks_interval_dict[i]
        conserverd_blocks_interval_dict_4.update(new_moved_dict_4)
        nonconserverd_blocks_interval_dict.pop(i)
    tmp4 = []
    for i in sorted(conserverd_blocks_interval_dict_4.keys()):
        tmp4.append(conserverd_blocks_interval_dict_4[i])
    result4 = pd.concat(tmp4)
    # print(result4)
    result4.to_excel(writer, sheet_name='4', header=None, index=False)


    #---------conserverd5----------------
    conserverd_df_5 = pd.read_excel(conserverd_group, sheet_name=4, names=['ROC_id','ROC_chr','ROC_group','ROC_gene_position','Sb_id',
                                                                               'Sb_chr','Sb_gene_position'],header=None)
    #----------拆分为每一块---------------
    conserverd_split_blocks_row_5 = []
    for index,row in conserverd_df_5.iterrows():
        if conserverd_df_5.loc[index]['ROC_id'] == '###':
            conserverd_split_blocks_row_5.append(index)
    conserverd_split_blocks_dict_5 = {}
    for i in range(0,len(conserverd_split_blocks_row_5)):
        if i != len(conserverd_split_blocks_row_5)-1:
            conserverd_split_blocks_dict_5[i] = conserverd_df_5.iloc[conserverd_split_blocks_row_5[i]:conserverd_split_blocks_row_5[i+1],:]
        else:
            conserverd_split_blocks_dict_5[i] = conserverd_df_5.iloc[conserverd_split_blocks_row_5[i]:conserverd_df_5.shape[0],:]
    # print(conserverd_split_blocks_dict)
    #-------------将保守的group按照max、min值构建区间-------------
    conserverd_blocks_interval_dict_5 = {}
    for i in conserverd_split_blocks_dict_5.keys():
        block_start = conserverd_split_blocks_dict_5[i]['Sb_gene_position'].min()
        block_end   = conserverd_split_blocks_dict_5[i]['Sb_gene_position'].max()
        conserverd_blocks_interval_dict_5[Interval(block_start,block_end)] = conserverd_split_blocks_dict_5[i]
    # print(conserverd_blocks_interval_dict_5)

    #----------------区间判断----------------------------
    conserverd_blocks_interval_intertree_5 = IntervalTree(conserverd_blocks_interval_dict_5)
    not_overlaps_inter_5 = []
    for inter in nonconserverd_blocks_interval_dict.keys():
        if not conserverd_blocks_interval_intertree_5.overlaps(inter):
            not_overlaps_inter_5.append(inter)
    # print(not_overlaps_inter_5)
    ##----把需要转移的block进行转移------------------------
    new_moved_dict_5 = {}
    for i in not_overlaps_inter_5:
        new_moved_dict_5[i] = nonconserverd_blocks_interval_dict[i]
        conserverd_blocks_interval_dict_5.update(new_moved_dict_5)
        nonconserverd_blocks_interval_dict.pop(i)
    tmp5 = []
    for i in sorted(conserverd_blocks_interval_dict_5.keys()):
        tmp5.append(conserverd_blocks_interval_dict_5[i])
    result5 = pd.concat(tmp5)
    result5.to_excel(writer, sheet_name='5', header=None, index=False)



    #---------conserverd6----------------
    conserverd_df_6 = pd.read_excel(conserverd_group, sheet_name=5, names=['ROC_id','ROC_chr','ROC_group','ROC_gene_position','Sb_id',
                                                                               'Sb_chr','Sb_gene_position'],header=None)
    #----------拆分为每一块---------------
    conserverd_split_blocks_row_6 = []
    for index,row in conserverd_df_6.iterrows():
        if conserverd_df_6.loc[index]['ROC_id'] == '###':
            conserverd_split_blocks_row_6.append(index)
    conserverd_split_blocks_dict_6 = {}
    for i in range(0,len(conserverd_split_blocks_row_6)):
        if i != len(conserverd_split_blocks_row_6)-1:
            conserverd_split_blocks_dict_6[i] = conserverd_df_6.iloc[conserverd_split_blocks_row_6[i]:conserverd_split_blocks_row_6[i+1],:]
        else:
            conserverd_split_blocks_dict_6[i] = conserverd_df_6.iloc[conserverd_split_blocks_row_6[i]:conserverd_df_6.shape[0],:]
    # print(conserverd_split_blocks_dict)
    #-------------将保守的group按照max、min值构建区间-------------
    conserverd_blocks_interval_dict_6 = {}
    for i in conserverd_split_blocks_dict_6.keys():
        block_start = conserverd_split_blocks_dict_6[i]['Sb_gene_position'].min()
        block_end   = conserverd_split_blocks_dict_6[i]['Sb_gene_position'].max()
        conserverd_blocks_interval_dict_6[Interval(block_start,block_end)] = conserverd_split_blocks_dict_6[i]
    # print(conserverd_blocks_interval_dict_6)

    #----------------区间判断----------------------------
    conserverd_blocks_interval_intertree_6 = IntervalTree(conserverd_blocks_interval_dict_6)
    not_overlaps_inter_6 = []
    for inter in nonconserverd_blocks_interval_dict.keys():
        if not conserverd_blocks_interval_intertree_6.overlaps(inter):
            not_overlaps_inter_6.append(inter)
    # print(not_overlaps_inter_6)
    ##----把需要转移的block进行转移------------------------
    new_moved_dict_6 = {}
    for i in not_overlaps_inter_6:
        new_moved_dict_6[i] = nonconserverd_blocks_interval_dict[i]
        conserverd_blocks_interval_dict_6.update(new_moved_dict_6)
        nonconserverd_blocks_interval_dict.pop(i)
    tmp6 = []
    for i in sorted(conserverd_blocks_interval_dict_6.keys()):
        tmp6.append(conserverd_blocks_interval_dict_6[i])
    result6 = pd.concat(tmp6)
    result6.to_excel(writer, sheet_name='6', header=None, index=False)



    #---------conserverd7----------------
    conserverd_df_7 = pd.read_excel(conserverd_group, sheet_name=6, names=['ROC_id','ROC_chr','ROC_group','ROC_gene_position','Sb_id',
                                                                               'Sb_chr','Sb_gene_position'],header=None)
    #----------拆分为每一块---------------
    conserverd_split_blocks_row_7 = []
    for index,row in conserverd_df_7.iterrows():
        if conserverd_df_7.loc[index]['ROC_id'] == '###':
            conserverd_split_blocks_row_7.append(index)
    conserverd_split_blocks_dict_7 = {}
    for i in range(0,len(conserverd_split_blocks_row_7)):
        if i != len(conserverd_split_blocks_row_7)-1:
            conserverd_split_blocks_dict_7[i] = conserverd_df_7.iloc[conserverd_split_blocks_row_7[i]:conserverd_split_blocks_row_7[i+1],:]
        else:
            conserverd_split_blocks_dict_7[i] = conserverd_df_7.iloc[conserverd_split_blocks_row_7[i]:conserverd_df_7.shape[0],:]
    # print(conserverd_split_blocks_dict)
    #-------------将保守的group按照max、min值构建区间-------------
    conserverd_blocks_interval_dict_7 = {}
    for i in conserverd_split_blocks_dict_7.keys():
        block_start = conserverd_split_blocks_dict_7[i]['Sb_gene_position'].min()
        block_end   = conserverd_split_blocks_dict_7[i]['Sb_gene_position'].max()
        conserverd_blocks_interval_dict_7[Interval(block_start,block_end)] = conserverd_split_blocks_dict_7[i]
    # print(conserverd_blocks_interval_dict_7)

    #----------------区间判断----------------------------
    conserverd_blocks_interval_intertree_7 = IntervalTree(conserverd_blocks_interval_dict_7)
    not_overlaps_inter_7 = []
    for inter in nonconserverd_blocks_interval_dict.keys():
        if not conserverd_blocks_interval_intertree_7.overlaps(inter):
            not_overlaps_inter_7.append(inter)
    # print(not_overlaps_inter_7)
    ##----把需要转移的block进行转移------------------------
    new_moved_dict_7 = {}
    for i in not_overlaps_inter_7:
        new_moved_dict_7[i] = nonconserverd_blocks_interval_dict[i]
        conserverd_blocks_interval_dict_7.update(new_moved_dict_7)
        nonconserverd_blocks_interval_dict.pop(i)
    tmp7 = []
    for i in sorted(conserverd_blocks_interval_dict_7.keys()):
        tmp7.append(conserverd_blocks_interval_dict_7[i])
    result7 = pd.concat(tmp7)
    result7.to_excel(writer, sheet_name='7', header=None, index=False)


    #---------conserverd8----------------
    conserverd_df_8 = pd.read_excel(conserverd_group, sheet_name=7, names=['ROC_id','ROC_chr','ROC_group','ROC_gene_position','Sb_id',
                                                                               'Sb_chr','Sb_gene_position'],header=None)
    #----------拆分为每一块---------------
    conserverd_split_blocks_row_8 = []
    for index,row in conserverd_df_8.iterrows():
        if conserverd_df_8.loc[index]['ROC_id'] == '###':
            conserverd_split_blocks_row_8.append(index)
    conserverd_split_blocks_dict_8 = {}
    for i in range(0,len(conserverd_split_blocks_row_8)):
        if i != len(conserverd_split_blocks_row_8)-1:
            conserverd_split_blocks_dict_8[i] = conserverd_df_8.iloc[conserverd_split_blocks_row_8[i]:conserverd_split_blocks_row_8[i+1],:]
        else:
            conserverd_split_blocks_dict_8[i] = conserverd_df_8.iloc[conserverd_split_blocks_row_8[i]:conserverd_df_8.shape[0],:]
    # print(conserverd_split_blocks_dict)
    #-------------将保守的group按照max、min值构建区间-------------
    conserverd_blocks_interval_dict_8 = {}
    for i in conserverd_split_blocks_dict_8.keys():
        block_start = conserverd_split_blocks_dict_8[i]['Sb_gene_position'].min()
        block_end   = conserverd_split_blocks_dict_8[i]['Sb_gene_position'].max()
        conserverd_blocks_interval_dict_8[Interval(block_start,block_end)] = conserverd_split_blocks_dict_8[i]
    # print(conserverd_blocks_interval_dict_8)

    #----------------区间判断----------------------------
    conserverd_blocks_interval_intertree_8 = IntervalTree(conserverd_blocks_interval_dict_8)
    not_overlaps_inter_8 = []
    for inter in nonconserverd_blocks_interval_dict.keys():
        if not conserverd_blocks_interval_intertree_8.overlaps(inter):
            not_overlaps_inter_8.append(inter)
    # print(not_overlaps_inter_8)
    ##----把需要转移的block进行转移------------------------
    new_moved_dict_8 = {}
    for i in not_overlaps_inter_8:
        new_moved_dict_8[i] = nonconserverd_blocks_interval_dict[i]
        conserverd_blocks_interval_dict_8.update(new_moved_dict_8)
        nonconserverd_blocks_interval_dict.pop(i)
    tmp8 = []
    for i in sorted(conserverd_blocks_interval_dict_8.keys()):
        tmp8.append(conserverd_blocks_interval_dict_8[i])
    result8 = pd.concat(tmp8)
    result8.to_excel(writer, sheet_name='8', header=None, index=False)

    #---------conserverd9----------------
    conserverd_df_9 = pd.read_excel(conserverd_group, sheet_name=8, names=['ROC_id','ROC_chr','ROC_group','ROC_gene_position','Sb_id',
                                                                               'Sb_chr','Sb_gene_position'],header=None)
    #----------拆分为每一块---------------
    conserverd_split_blocks_row_9 = []
    for index,row in conserverd_df_9.iterrows():
        if conserverd_df_9.loc[index]['ROC_id'] == '###':
            conserverd_split_blocks_row_9.append(index)
    conserverd_split_blocks_dict_9 = {}
    for i in range(0,len(conserverd_split_blocks_row_9)):
        if i != len(conserverd_split_blocks_row_9)-1:
            conserverd_split_blocks_dict_9[i] = conserverd_df_9.iloc[conserverd_split_blocks_row_9[i]:conserverd_split_blocks_row_9[i+1],:]
        else:
            conserverd_split_blocks_dict_9[i] = conserverd_df_9.iloc[conserverd_split_blocks_row_9[i]:conserverd_df_9.shape[0],:]
    # print(conserverd_split_blocks_dict)
    #-------------将保守的group按照max、min值构建区间-------------
    conserverd_blocks_interval_dict_9 = {}
    for i in conserverd_split_blocks_dict_9.keys():
        block_start = conserverd_split_blocks_dict_9[i]['Sb_gene_position'].min()
        block_end   = conserverd_split_blocks_dict_9[i]['Sb_gene_position'].max()
        conserverd_blocks_interval_dict_9[Interval(block_start,block_end)] = conserverd_split_blocks_dict_9[i]
    # print(conserverd_blocks_interval_dict_9)

    #----------------区间判断----------------------------
    conserverd_blocks_interval_intertree_9 = IntervalTree(conserverd_blocks_interval_dict_9)
    not_overlaps_inter_9 = []
    for inter in nonconserverd_blocks_interval_dict.keys():
        if not conserverd_blocks_interval_intertree_9.overlaps(inter):
            not_overlaps_inter_9.append(inter)
    ##----把需要转移的block进行转移------------------------
    new_moved_dict_9 = {}
    for i in not_overlaps_inter_9:
        new_moved_dict_9[i] = nonconserverd_blocks_interval_dict[i]
        conserverd_blocks_interval_dict_9.update(new_moved_dict_9)
        nonconserverd_blocks_interval_dict.pop(i)
    tmp9 = []
    for i in sorted(conserverd_blocks_interval_dict_9.keys()):
        tmp9.append(conserverd_blocks_interval_dict_9[i])
    result9 = pd.concat(tmp9)
    result9.to_excel(writer, sheet_name='9', header=None, index=False)


    #---------conserverd10----------------
    conserverd_df_10 = pd.read_excel(conserverd_group, sheet_name=9, names=['ROC_id','ROC_chr','ROC_group','ROC_gene_position','Sb_id',
                                                                               'Sb_chr','Sb_gene_position'],header=None)
    #----------拆分为每一块---------------
    conserverd_split_blocks_row_10 = []
    for index,row in conserverd_df_10.iterrows():
        if conserverd_df_10.loc[index]['ROC_id'] == '###':
            conserverd_split_blocks_row_10.append(index)
    conserverd_split_blocks_dict_10 = {}
    for i in range(0,len(conserverd_split_blocks_row_10)):
        if i != len(conserverd_split_blocks_row_10)-1:
            conserverd_split_blocks_dict_10[i] = conserverd_df_10.iloc[conserverd_split_blocks_row_10[i]:conserverd_split_blocks_row_10[i+1],:]
        else:
            conserverd_split_blocks_dict_10[i] = conserverd_df_10.iloc[conserverd_split_blocks_row_10[i]:conserverd_df_10.shape[0],:]
    # print(conserverd_split_blocks_dict)
    #-------------将保守的group按照max、min值构建区间-------------
    conserverd_blocks_interval_dict_10 = {}
    for i in conserverd_split_blocks_dict_10.keys():
        block_start = conserverd_split_blocks_dict_10[i]['Sb_gene_position'].min()
        block_end   = conserverd_split_blocks_dict_10[i]['Sb_gene_position'].max()
        conserverd_blocks_interval_dict_10[Interval(block_start,block_end)] = conserverd_split_blocks_dict_10[i]
    # print(conserverd_blocks_interval_dict_10)

    #----------------区间判断----------------------------
    conserverd_blocks_interval_intertree_10 = IntervalTree(conserverd_blocks_interval_dict_10)
    not_overlaps_inter_10 = []
    for inter in nonconserverd_blocks_interval_dict.keys():
        if not conserverd_blocks_interval_intertree_10.overlaps(inter):
            not_overlaps_inter_10.append(inter)
    ##----把需要转移的block进行转移------------------------
    new_moved_dict_10 = {}
    for i in not_overlaps_inter_10:
        new_moved_dict_10[i] = nonconserverd_blocks_interval_dict[i]
        conserverd_blocks_interval_dict_10.update(new_moved_dict_10)
        nonconserverd_blocks_interval_dict.pop(i)
    tmp10 = []
    for i in sorted(conserverd_blocks_interval_dict_10.keys()):
        tmp10.append(conserverd_blocks_interval_dict_10[i])
    result10 = pd.concat(tmp10)
    result10.to_excel(writer, sheet_name='10', header=None, index=False)
    writer.save()

    remaining_non_group_list = []
    for l in nonconserverd_blocks_interval_dict.keys():
        # print(l)
        # print(nonconserverd_blocks_interval_dict[l])
        remaining_non_group_list.append(nonconserverd_blocks_interval_dict[l])
    if remaining_non_group_list:
        remaining_non_group = pd.concat(remaining_non_group_list)
        remaining_non_group.to_excel('remain' + non_group_id + '.xlsx', header=False, index=False)


def get_run_script(xlsx1,xlsx2,xlsx3,xlsx4,xlsx5,Output_script):
    with open(Output_script,'w') as f:
        xlsx1_20 = xlrd.open_workbook(xlsx1)
        xlsx1_20_index = len(xlsx1_20.sheets())
        f.write('#'+'xlsx1-20 total sheets num : {}'.format(str(xlsx1_20_index))+'\n')
        f.write('python Roc22_chrCorrect.py -c conserved_group.xlsx -n nonconserved_group_1-20.xlsx -s 0 -o 0.xlsx'+'\n')
        for i in range(0,xlsx1_20_index-1):
            f.write('python Roc22_chrCorrect.py -c {}.xlsx -n nonconserved_group_1-20.xlsx -s {} -o {}.xlsx'.format(i,i+1,i+1)+'\n')

        xlsx21_40 = xlrd.open_workbook(xlsx2)
        xlsx21_40_index = len(xlsx21_40.sheets())
        f.write('#' + 'xlsx21-40 total sheets num : {}'.format(str(xlsx21_40_index))+'\n')
        for i in range(0+xlsx1_20_index-1,xlsx21_40_index+xlsx1_20_index-1):
            f.write('python Roc22_chrCorrect.py -c {}.xlsx -n nonconserved_group_21-40.xlsx -s {} -o {}.xlsx'.format(i,i-xlsx1_20_index+1,i+1)+'\n')


        xlsx41_60 = xlrd.open_workbook(xlsx3)
        xlsx41_60_index = len(xlsx41_60.sheets())
        f.write('#' + 'xlsx41-60 total sheets num : {}'.format(str(xlsx41_60_index))+'\n')
        for i in range(0 + xlsx1_20_index+xlsx21_40_index-1 , xlsx41_60_index + xlsx1_20_index +xlsx21_40_index-1):
            f.write('python Roc22_chrCorrect.py -c {}.xlsx -n nonconserved_group_41-60.xlsx -s {} -o {}.xlsx'.format(i, i -xlsx1_20_index-xlsx21_40_index+1,i + 1)+'\n')

        xlsx61_80 = xlrd.open_workbook(xlsx4)
        xlsx61_80_index = len(xlsx61_80.sheets())
        f.write('#' + 'xlsx61-80 total sheets num : {}'.format(str(xlsx61_80_index))+'\n')
        for i in range(0 + xlsx1_20_index + xlsx21_40_index + xlsx41_60_index-1, xlsx61_80_index + xlsx1_20_index + xlsx21_40_index+xlsx41_60_index-1):
            f.write('python Roc22_chrCorrect.py -c {}.xlsx -n nonconserved_group_61-80.xlsx -s {} -o {}.xlsx'.format(i, i -xlsx1_20_index - xlsx21_40_index - xlsx41_60_index+1, i + 1)+'\n')

        xlsx81_end = xlrd.open_workbook(xlsx5)
        xlsx81_end_index = len(xlsx81_end.sheets())
        f.write('#' + 'xlsx81-end total sheets num : {}'.format(str(xlsx81_end_index))+'\n')
        for i in range(0 + xlsx1_20_index + xlsx21_40_index + xlsx41_60_index + xlsx61_80_index-1,
                       xlsx81_end_index + xlsx1_20_index + xlsx21_40_index + xlsx41_60_index + xlsx61_80_index-1):
            f.write('python Roc22_chrCorrect.py -c {}.xlsx -n nonconserved_group_81-end.xlsx -s {} -o {}.xlsx'.format(i, i - xlsx1_20_index - xlsx21_40_index - xlsx41_60_index - xlsx61_80_index+1,i + 1)+'\n')

    f.close()


def main(args):
    if args.group and args.querylist:
        group_file = args.group
        query_list = args.querylist
        Separate_groupAndSort_blocks(group_file,query_list)
    if args.conserverd and args.nonconserverd:
        conserverd_group = args.conserverd
        nonconserverd_group = args.nonconserverd
        nonconserverd_sheet = args.sheetid
        output_file = args.output
        correct_chr(conserverd_group, nonconserverd_group, nonconserverd_sheet, output_file)
    if args.xlsx1 and args.xlsx2:
        xlsx1 = args.xlsx1
        xlsx2 = args.xlsx2
        xlsx3 = args.xlsx3
        xlsx4 = args.xlsx4
        xlsx5 = args.xlsx5
        Output_script = args.Outputscript
        get_run_script(xlsx1,xlsx2,xlsx3,xlsx4,xlsx5,Output_script)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='Roc22_chrCorrect.py',
        formatter_class=argparse.RawTextHelpFormatter,
        description='''This script is used for ROC22 genome chromosome correction work''',
        epilog= 'author:\t{0}\nmail:\t{1}\nversion:\t{2}'.format(__author__,__mail__,__version__))
    parser.add_argument('-g', '--group',  help='Input the uncorrected chromosome group file')
    parser.add_argument('-q', '--querylist', nargs='+', help='Input the id of the group to keep (eg:1,11,21,31...)')

    parser.add_argument('-c', '--conserverd',  help='Input conservative block file(.xlsx)')
    parser.add_argument('-n', '--nonconserverd', help='Input non-conservative block file(.xlsx)')
    parser.add_argument('-s', '--sheetid', type=int,  help='Input non-conservative file sheetID')
    parser.add_argument('-o', '--output', help='Output result file(.xlsx)')

    parser.add_argument('-x1', '--xlsx1', help='Input nonconserverd_group_1-20(.xlsx)')
    parser.add_argument('-x2', '--xlsx2', help='Input nonconserverd_group_21-40(.xlsx)')
    parser.add_argument('-x3', '--xlsx3', help='Input nonconserverd_group_41-60(.xlsx)')
    parser.add_argument('-x4', '--xlsx4', help='Input nonconserverd_group_61-80(.xlsx)')
    parser.add_argument('-x5', '--xlsx5', help='Input nonconserverd_group_81-end(.xlsx)')
    parser.add_argument('-O', '--Outputscript', help='Output shell script(.sh)')

    args = parser.parse_args()
    main(args)
