#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/3/12

'''
将保守的group与非保守的group的blocks全部读入到字典中，字典格式：groupID-indexID:con_blocks/non_blocks
'''

import pandas as pd
import re
from path import Path
from intervaltree import Interval, IntervalTree
from collections import OrderedDict


def Store_all_conblocks(con_df):
    # 将所有保守group的所有blocks存入字典，字典格式为：groupID:all_con_block_df    [ex: 7:all_group7_block_df]
    group2conblocksDict = {}
    for i in range(0,10):
        try:
            con = pd.read_excel(con_df,sheet_name=i,header=None,names=['ROC_id','ROC_chr','ROC_group','ROC_gene_position','Sb_id','Sb_chr','Sb_gene_position'])
            con_group_id = re.findall('[0-9]+', con.iloc[1, :]['ROC_group'])[0]
            group2conblocksDict[con_group_id] = con
        except:
            continue
    return group2conblocksDict


def Store_all_conblocks2(con_df):
    # 将所有保守group的blocks存入字典，字典格式为：groupID-indexID:block_df    [ex: 7-0:con_block_df]
    GroupIndex2blocksDict = {}
    for i in range(0,10):
        try:
            con = pd.read_excel(con_df,sheet_name=i,header=None,names=['ROC_id','ROC_chr','ROC_group','ROC_gene_position','Sb_id','Sb_chr','Sb_gene_position'])
            con_group_id = re.findall('[0-9]+', con.iloc[1, :]['ROC_group'])[0]
            con_split_blocks_row = []
            for index,row in con.iterrows():
                if con.loc[index]['ROC_id'] == '###':
                    con_split_blocks_row.append(index)

            for i in range(0,len(con_split_blocks_row)):
                if i != len(con_split_blocks_row)-1:
                    GroupIndex2blocksDict[con_group_id+'-'+str(i)] = con.iloc[con_split_blocks_row[i]:con_split_blocks_row[i+1],:]
                else:
                    GroupIndex2blocksDict[con_group_id+'-'+str(i)] = con.iloc[con_split_blocks_row[i]:con.shape[0],:]
        except:
            continue
    return GroupIndex2blocksDict


def Store_all_nonblocks(file_path):
    # 将所有保守group的blocks存入字典，字典格式为：groupID-indexID:block_df ex: 7-0:non_block_df
    files = [file for file in Path(file_path).files() if 'nonconserverd_group' in file.basename()]
    group_index2nonblocksDict = {}
    for file in files:
        for i in range(0, 20):
            try:
                non = pd.read_excel(file, sheet_name=i, header=None,names=['ROC_id', 'ROC_chr', 'ROC_group', 'ROC_gene_position', 'Sb_id', 'Sb_chr','Sb_gene_position'])
                # print(non)
                non_group_id = re.findall('[0-9]+', non.iloc[1, :]['ROC_group'])[0]
                # print(non_group_id)
                con_split_blocks_row = []

                for index, row in non.iterrows():
                    if non.loc[index]['ROC_id'] == '###':
                        con_split_blocks_row.append(index)

                for i in range(0, len(con_split_blocks_row)):
                    if i != len(con_split_blocks_row) - 1:
                        group_index2nonblocksDict[non_group_id + '-' + str(i)] = non.iloc[con_split_blocks_row[i]:con_split_blocks_row[i + 1], :]
                    else:
                        group_index2nonblocksDict[non_group_id + '-' + str(i)] = non.iloc[con_split_blocks_row[i]:non.shape[0],:]
            except:
                continue
    # for key in group_index2nonblocksDict.keys():
    #     print(key)
    #     print(group_index2nonblocksDict[key])
    # print(len(group_index2nonblocksDict.keys())) #450
    return group_index2nonblocksDict


def Store_conblocks2Interval(ConDict2):
    ## 将所有保守group的blocks存入有序字典，与输入字典格式一致，为：groupid-index:block_df   ex: 6-0:con_block_df
    ConSortblockDict = OrderedDict()
    for i in sorted(ConDict2.items(), key=lambda x: len(x[1]), reverse=True):
        ConSortblockDict[i[0]] = i[1]

    # 使用多键值，将Con_blocks按照max、min值构建区间保存为字典。字典格式 (Interval(min,max),group-index):non_block
    # ex:(Interval(45062.0, 943806.0), '6-0'):non_block_df
    con_blocks_interval_dict = OrderedDict()
    for key in ConSortblockDict.keys():
        block_start = ConSortblockDict[key]['Sb_gene_position'].min()
        block_end = ConSortblockDict[key]['Sb_gene_position'].max()
        con_blocks_interval_dict[(Interval(block_start, block_end), key)] = ConSortblockDict[key]
    return con_blocks_interval_dict


def Store_nonblocks2Interval(NonblockDict):
    '''
    :param NonblockDict: groupid-index:block_df   ex: 29-4:non_block_df
    :return: (Interval(min,max),group-index):non_block   ex:(Interval(58145965.0, 59310987.0), '29-4'):non_block_df
    '''
    ## 将所有非保守group的blocks存入有序字典，与输入字典格式一致，为：groupid-index:block_df   ex: 29-4:non_block_df
    NonSortblockDict = OrderedDict()
    for i in sorted(NonblockDict.items(),key=lambda x:len(x[1]),reverse=True):
        NonSortblockDict[i[0]] = i[1]

    # 使用多键值，将non_blocks按照max、min值构建区间保存为字典。字典格式 (Interval(min,max),group-index):non_block
    # ex:(Interval(58145965.0, 59310987.0), '29-4'):non_block_df
    non_blocks_interval_dict = OrderedDict()
    for key in NonSortblockDict.keys():
        block_start = NonSortblockDict[key]['Sb_gene_position'].min()
        block_end   = NonSortblockDict[key]['Sb_gene_position'].max()
        non_blocks_interval_dict[(Interval(block_start,block_end),key)] = NonSortblockDict[key]

    return non_blocks_interval_dict



def Correct_group5(ConinterDict,NoninterDict):

    #测试原输入文件中有没有下面需要删除或添加的block
    # for i in NoninterDict.keys():
    #     # print(i)
    #     if str(i) == "(Interval(58145965.0, 59310987.0), '29-4')":
    #         print(i)
    # for i in ConinterDict.keys():
    #     if str(i) == "(Interval(57947569.0, 58132767.0), '70-10')":
    #         print(i)
    # print('-'*20)

    print(len(ConinterDict.keys()))
    print(len(NoninterDict.keys()))

    #去除排好顺序后NoninterDict中的第一个block，即最长的block
    firstblock = list(NoninterDict.items())[0]
    # print(firstblock[0])
    # print(firstblock[1])


    #将读入的ConinterDict按照Con的group分成10个不同的dict
    All_sub_ConinterDict = {}
    ConGroup = [1,11,12,19,22,24,59,73]

    for Conblock in ConinterDict.keys():
        # print(Conblock)
        # print(ConinterDict[Conblock])
        for num in ConGroup:
            if Conblock[1].split('-')[0] == str(num):
                if Conblock[1].split('-')[0] not in All_sub_ConinterDict.keys():
                    All_sub_ConinterDict[Conblock[1].split('-')[0]] = {Conblock:ConinterDict[Conblock]}
                else:
                    All_sub_ConinterDict[Conblock[1].split('-')[0]].update({Conblock:ConinterDict[Conblock]})

    ConNot_overlap_Conblock = OrderedDict()  # 先检查一下这个Nonblock与Conblock是否有overlap，把没有overlap的输出到本字典中

    Notoverlap_closeblock = OrderedDict()  # 判断这个Nonblock与哪个group没有overlap，然后针对这几个没有overlap的group中选出最接近的block
    for i in All_sub_ConinterDict.keys(): #i=6,7,8,13,...
        list2tree = []
        # print(i)  #7
        # print(All_sub_ConinterDict[i].keys())   #dict_keys([(Interval(6161098.0, 47279050.0), '7-5'), (Interval(9358715.0, 44029596.0), '7-4'),...
        for j in All_sub_ConinterDict[i].keys():
            list2tree.append(j[0])
        # print(list2tree)
        CongroupTree = IntervalTree(list2tree)
        # print(CongroupTree)

        if not CongroupTree.overlap(firstblock[0][0]):
            # ConNot_overlap_Conblock[All_sub_ConinterDict[i]] = ConinterDict[i]
            # print(i)
            Lscore = OrderedDict()
            for key in All_sub_ConinterDict[i]:
                # print(key)
                if key[0].end <= firstblock[0][0].begin:
                    Lscore[key] = firstblock[0][0].begin - key[0].end
                elif firstblock[0][0].end <= key[0].begin:
                    Lscore[key] = key[0].begin - firstblock[0][0].end
                # break
            # for i in Lscore:
            #     print(i)
            #     print(Lscore[i])
            #     print('*'*50)
            group_closestNonblock = list(sorted(Lscore.items(), key=lambda x: x[1], reverse=False))[0]
            # print(group_closestNonblock)  #输出没有overlap的group中与nonblock最近的一组block
            Notoverlap_closeblock[i] = group_closestNonblock
        else:
            continue
        # print('-' * 200)

    # for i in Notoverlap_closeblock.keys():
    #     print(i)
    #     print(Notoverlap_closeblock[i])
    #     break
    best_closestblock = list(sorted(Notoverlap_closeblock.items(), key=lambda x: x[1][1], reverse=False))[0]
    # print(best_closestblock)
    print('{} need move to {}'.format(firstblock[0][1], best_closestblock[1][0][1]))

    # print(firstblock[0])
    new_name = (firstblock[0][0],best_closestblock[1][0][1].split('-')[0]+'-'+'new') #对将要移动的nonblock进行重命名，命名为移动到目的地group的 groupID-new
    # print(new_name)
    ConinterDict.update({new_name:NoninterDict[firstblock[0]]})
    NoninterDict.pop(firstblock[0])

    print(len(ConinterDict.keys()))
    print(len(NoninterDict.keys()))

    return ConinterDict,NoninterDict




if __name__ == '__main__':
    # conDict = Store_all_conblocks(r'E:\Badila\synteny_analysis\C03\parsed1\conserverd_group.xlsx')
    conDict2 = Store_all_conblocks2(r'E:\Badila\synteny_analysis\C03\conserverd_group.xlsx')
    nonDict = Store_all_nonblocks(r'E:\Badila\synteny_analysis\C03')
    nonIntervalDict = Store_nonblocks2Interval(nonDict)
    conIntervalDict = Store_conblocks2Interval(conDict2)
    print(len(conDict2.keys()))
    print(len(nonDict.keys()))
    print(len(nonIntervalDict.keys()))
    print(len(conIntervalDict.keys()))
