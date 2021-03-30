#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/3/12

__author__ = 'Haoran Pan'
__mail__ = 'haoran_pan@qq.com'
__date__ = '20210330'
__version__ = 'v3.0'

import pandas as pd
import CPCS1
import re
import sys
import argparse


def parse_ConGroups(group_df,conserverd_group_num,conserverd_group_concat):
    writer_conserverd = pd.ExcelWriter('conserverd_group.xlsx')
    for num in conserverd_group_num:
        conserverd_group_individual_concat = []
        for j in conserverd_group_concat:
            group_id = re.findall('[0-9]+',j.iloc[1,:]['ROC_group'])[0]
            if int(group_id) == int(num):
                conserverd_group_individual_concat.append(j)
        conserverd_group_individual=pd.concat(conserverd_group_individual_concat)
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


def parse_NonGroups(group_df,conserverd_group_num,nonconserverd_group_concat):
    split_groups_num = split_num()


    for split_group in split_groups_num[:-1]:
        writer_nonconserverd = pd.ExcelWriter('nonconserverd_group_'+str(split_group[0])+'_'+str(split_group[-1])+'.xlsx')
        for num in split_group:
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

                # --------------------------------------------------------------
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


def main(args):
    sort_tab = args.input
    con_list = args.conserved

    group_df = CPCS1.remove_redundant_comment_line(sort_tab)
    ConGroup_num, non, con = CPCS1.SeparateClassifiBlocks(group_df,con_list)
    parse_ConGroups(group_df, ConGroup_num, con)
    split_num()
    parse_NonGroups(group_df, ConGroup_num, non)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Delete the redundant "#" between segments in each group, split each segment according to "#", and divide the conservative and non-conservative ones into a dataframe respectively.\nThe 10 conservative groups and the remaining non-conservative groups are sorted separately and output as xlsx files''',
        usage="python {} -i <sort.xlsx> -c conserved_group_num ".format(sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__, __mail__, __date__,__version__))
    parser.add_argument('-i','--input',required=True,help='Input sorted xlsx')
    parser.add_argument('-c','--conserved',nargs='+',required=True,help='Input conserved group number')
    args = parser.parse_args()
    main(args)
