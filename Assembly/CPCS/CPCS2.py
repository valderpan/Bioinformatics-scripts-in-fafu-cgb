#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/3/12

'''
将10组保守的group与剩余非保守的group分别进行整合排序后输出为xlsx文件
'''

import pandas as pd
import CPCS1
import re

pd.set_option('display.max_columns', 10)

def parse_ConGroups(group_df,conserverd_group_num,conserverd_group_concat):
    '''
    排序并输出10组保守的group
    :param conserverd_file:将10组gruop保存在一个excel文件中
    :return:
    '''
    print(group_df)
    writer_conserverd = pd.ExcelWriter('conserverd_group.xlsx')
    for num in conserverd_group_num: #[6,7,8,13,31,49,62,70,73,103]
        # print(num)
        #--------将每个保守的group block整合后分别输出-----------------------
        conserverd_group_individual_concat = []
        for j in conserverd_group_concat:
            group_id = re.findall('[0-9]+',j.iloc[1,:]['ROC_group'])[0]
            # print(group_id)
            # break
            if int(group_id) == int(num):
                # print(j)
                conserverd_group_individual_concat.append(j)
        # print(conserverd_group_individual_concat)
        conserverd_group_individual=pd.concat(conserverd_group_individual_concat)
        #--------------------------------------------------------------
        #-------------对每个group进行'###'分割为blocks，排序后整合输出到xlsx中---------------------
        sub_split_blocks_row = []
        for index, row in conserverd_group_individual.iterrows():
            if group_df.loc[index]['ROC_id'] == '###':
                sub_split_blocks_row.append(index)

        group_split_blocks_dict = {}
        for i in range(0, len(sub_split_blocks_row)):
            if i != len(sub_split_blocks_row) - 1:
                group_split_blocks_dict[i] = group_df.iloc[sub_split_blocks_row[i]:sub_split_blocks_row[i + 1], :]
            else:
                group_split_blocks_dict[i] = group_df.iloc[
                                             sub_split_blocks_row[i]:conserverd_group_individual.index[-1] + 1, :]
        #这里是按照每个blocks中高粱gene的位置均值来排序
        con2sort = {}
        for j in group_split_blocks_dict.keys():
            con2sort[int(group_split_blocks_dict[j]['Sb_gene_position'].mean())] = group_split_blocks_dict[j]

        sorted_blocks_concat = []
        for i in sorted(con2sort):
            sorted_blocks_concat.append(con2sort[i])
        sorted_blocks = pd.concat(sorted_blocks_concat)

        sorted_blocks.to_excel(writer_conserverd, sheet_name=str(num), header=None, index=False)
    writer_conserverd.save()


def split_num():
    all_groups_num = [i for i in range(1,80)]   #80个group就将range改为1,80 #TODO
    split_group_num = [all_groups_num[i:i+20] for i in range(0,80,20)]
    return split_group_num
    # print(split_group_num)
    # for i in split_group_num:
    #     print(len(i))

def parse_NonGroups(group_df,conserverd_group_num,nonconserverd_group_concat):
    split_groups_num = split_num()
    # print(split_groups_num)
    # group_df = remove_redundant_comment_line(group_file)
    # group_df = group_df.reset_index(drop=True)

    for split_group in split_groups_num[:-1]:
        writer_nonconserverd = pd.ExcelWriter('nonconserverd_group_'+str(split_group[0])+'_'+str(split_group[-1])+'.xlsx')
        for num in split_group:
            if str(num) not in conserverd_group_num:
                # --------将每个保守的group block整合后分别输出-----------------------
                nonconserverd_group_individual_concat = []
                for j in nonconserverd_group_concat:
                    group_id = re.findall('[0-9]+', j.iloc[1, :]['ROC_group'])[0]
                    if group_id == str(num):
                        nonconserverd_group_individual_concat.append(j)
                try:
                    nonconserverd_group_individual = pd.concat(nonconserverd_group_individual_concat)
                except:
                    continue

                # --------------------------------------------------------------
                # -------------对每个group进行'###'分割为blocks，排序后整合输出到xlsx中---------------------
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
                # 这里是按照每个blocks中高粱gene的位置均值来排序
                non2sort_1 = {}
                for j in group_split_blocks_dict.keys():
                    non2sort_1[int(group_split_blocks_dict[j]['Sb_gene_position'].mean())] = group_split_blocks_dict[j]

                sorted_blocks_concat = []
                for i in sorted(non2sort_1):
                    sorted_blocks_concat.append(non2sort_1[i])
                sorted_blocks = pd.concat(sorted_blocks_concat)
                sorted_blocks.to_excel(writer_nonconserverd, sheet_name=str(num), header=None, index=False)
        writer_nonconserverd.save()

    last_split_group =  split_groups_num[-1]
    last_writer = pd.ExcelWriter('nonconserverd_group_'+str(last_split_group[0])+'_'+'end'+'.xlsx')
    for num in last_split_group:
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
            non2sort_2 = {}
            for j in group_split_blocks_dict.keys():
                non2sort_2[int(group_split_blocks_dict[j]['Sb_gene_position'].mean())] = group_split_blocks_dict[j]

            sorted_blocks_concat = []
            for i in sorted(non2sort_2):
                sorted_blocks_concat.append(non2sort_2[i])
            sorted_blocks = pd.concat(sorted_blocks_concat)
            sorted_blocks.to_excel(last_writer, sheet_name=str(num), header=None, index=False)
    last_writer.save()



if __name__ == '__main__':
    group_df = CPCS1.remove_redundant_comment_line(r'E:\Badila\synteny_analysis\C03\parsed1\C03.parse1.xlsx')
    ConGroup_num,non,con = CPCS1.SeparateClassifiBlocks(group_df)
    parse_ConGroups(group_df,ConGroup_num,con)
    split_num()
    parse_NonGroups(group_df, ConGroup_num, non)
