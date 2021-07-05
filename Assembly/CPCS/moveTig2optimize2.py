#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/6/11

'''
%prog <Origin_groupfile_path>

此脚本用来解完麻子后，根据movestat表将groupXX.txt中的congtig在各个group之间进行移动，并最后生成新的groupXX.txt
'''

import pandas as pd
from path import Path
import re
import sys
sys.path.append(r'D:\pycharm\panhr\FAFU_CGB\Genomic_exercise')
import SetLog


def read_origin_group(group_path):
    '''
    读取初始的所有的group.txt中的tig信息合并为一个总表;将{groupID tigID:他们所在行的dataframe}存储到字典
    :param group_path:
    :return:
    1.allOrigingroup_df
    2.Dict(group tig:row) eg : group1	tig00000047_0:0  tig00000047_0       624  1635221  group1
    '''
    files = [file for file in Path(group_path).files() if file.endswith('.txt')]
    concat_df = []
    grouptig2df = {}


    for file in files:
        file_name = re.findall('(group[0-9]+)\.txt', file.basename())[0]
        df = pd.read_table(file, comment='#', names=['Contig', 'RECounts', 'Length'])
        df['Group'] = file_name

        concat_df.append(df)
        for index, row in df.iterrows():  # TODO
            grouptig2df['{}'.format(row['Contig'])] = df.loc[[index]]
    allOrigingroup_df = pd.concat(concat_df, ignore_index=True)


    if len(grouptig2df) == allOrigingroup_df.shape[0]:  # TODO   会有一个contig在不同的group中;最初的assembly结果是没有这种情况
        SetLog.info_out('There is no duplicate tig in allOrigingroup dataframe! | [check ok ~]')
    else:
        SetLog.error_out('AllOrigingroup dataframe contains duplicate tig!')
        SetLog.error_out('Allorigingroup tig number : {}'.format(allOrigingroup_df.shape[0]))
        SetLog.error_out('grouptig2df key number : {}'.format(len(grouptig2df)))
    return allOrigingroup_df, grouptig2df


def read_file(movestat_xlsx):
    SetLog.info_out('Read movestat xlsx')
    movestat_df = pd.read_excel(movestat_xlsx,names=['tigID','pos1','pos2'])
    # print(movestat_df)
    return movestat_df


def Reorganize_group(Allorigroup, movestat_df, output_path=0):
    '''
    对原有的group进行重排，重拍后的groupxx.txt = origin - remove + add
    :param Allorigroup:
    :param movestat_df:
    :param output_path:
    :return:
    '''

    for num in range(1, 81):
        if num != 81:
            Oridf_query_df = Allorigroup[Allorigroup['Group'] == 'group{}'.format(num)]

            Remove_query_df = movestat_df[movestat_df['pos1'] == 'group{}'.format(num)]

            Remove_query_ctg = Remove_query_df['tigID'].tolist()

            Add_query_df = movestat_df[movestat_df['pos2'] == 'group{}'.format(num)]

            Add_query_ctg = Add_query_df['tigID'].tolist()

            concat_list = []
            removed_query_df = Oridf_query_df[~Oridf_query_df['Contig'].isin(Remove_query_ctg)]

            needadd_query_df = Allorigroup[Allorigroup['Contig'].isin(Add_query_ctg)]

            result_query_df = pd.concat([removed_query_df,needadd_query_df],ignore_index=True)

            result_query_df = result_query_df.iloc[:,[0,1,2]]
            result_query_df.rename(columns={'Contig':'#Contig'},inplace=True)

            if Oridf_query_df.shape[0] - len(Remove_query_ctg) + len(Add_query_ctg) == result_query_df.shape[0]:
                SetLog.info_out('Group{} || The number of ctg moved is consistent with the number of ctg counted ~'.format(num))
                if output_path:
                    result_query_df.drop_duplicates('#Contig',keep='first',inplace = True)
                    result_query_df.to_csv('{}/group{}.txt'.format(output_path,num),sep='\t',header=True,index=False)

            else:
                SetLog.error_out('Group{} || The number of ctg moved is inconsistent with the number of ctg counted !!!'.format(num))
                SetLog.error_out('Group{} Origin ctg : {}'.format(num,Oridf_query_df.shape[0]))
                SetLog.error_out('Group{} need remove ctg : {}'.format(num,len(Remove_query_ctg)))
                SetLog.error_out('Group{} need add ctg : {}'.format(num,len(Add_query_ctg)))
                SetLog.error_out('Group{} result ctg should is {} but it is {}'.format(num,Oridf_query_df.shape[0] - len(Remove_query_ctg) + len(Add_query_ctg),result_query_df.shape[0]))
                for ctg in Remove_query_df['tigID'].tolist():
                    if ctg not in Oridf_query_df['Contig'].tolist():
                        SetLog.warning_out('Group{} {} in Remove but it not in Origin'.format(num,ctg))
                if output_path:
                    result_query_df.drop_duplicates('#Contig', keep='first',inplace = True)
                    result_query_df.to_csv('{}/group{}.txt'.format(output_path,num),sep='\t',header=True,index=False)



if __name__ == '__main__':
    allorigingroup_df,grouptig2df = read_origin_group(r'E:\Badila\synteny_analysis\Badila染色体第二次校正\矫正后\afterround1\changename')
    #allorigingroup_df,grouptig2df = read_origin_group(r'E:\Badila\synteny_analysis\Badila染色体第一次校正\矫正后\movetigstat\origin_group_test') #最初的 一次也没矫正过的
    allorigingroup_df.to_excel('allorigingroup_afterround1df.xlsx',header=False,index=False)
    # allorigingroup_df,grouptig2df = read_origin_group(r'E:\Badila\synteny_analysis\Badila染色体第一次校正\矫正后\movetigstat\firstcorrected_group')
    # print(allorigingroup_df)
    # needmove_df = need2movetig(r'E:\Badila\synteny_analysis\Badila染色体第二次校正\矫正后\ZG_correct_round2.tigmovestat.allchr.xlsx',10)
    movestat_df = read_file(r'E:\Badila\synteny_analysis\Badila染色体第二次校正\矫正后\movestat.round2.clean.xlsx')
    Reorganize_group(allorigingroup_df,movestat_df,r'E:\Badila\synteny_analysis\Badila染色体第二次校正\矫正后\afterround2\changeID')
    # Reorganize_group(allorigingroup_df,needmove_df,r'E:\Badila\synteny_analysis\Badila染色体第二次校正\矫正后\corrected_groups_round2')

    # allorigingroup_df, grouptig2df = read_origin_group(sys.argv[1])
    # print(allorigingroup_df)
    # needmove_df = need2movetig(sys.argv[2], sys.argv[3])
    # Reorganize_group(allorigingroup_df, needmove_df)
