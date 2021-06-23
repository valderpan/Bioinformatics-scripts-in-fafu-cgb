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
            # grouptig2df['{}\t{}'.format(file_name, row['Contig'])] = df.loc[[index]]
            grouptig2df['{}'.format(row['Contig'])] = df.loc[[index]]
    allOrigingroup_df = pd.concat(concat_df, ignore_index=True)
    # for key in grouptig2df.keys():
    #     print(key)
    #     print(grouptig2df[key])
    #     break

    if len(grouptig2df) == allOrigingroup_df.shape[0]:  # TODO   会有一个contig在不同的group中
        SetLog.info_out('There is no duplicate tig in allOrigingroup dataframe! | [check ok ~]')
    else:
        SetLog.error_out('AllOrigingroup dataframe contains duplicate tig!')
        SetLog.error_out('Allorigingroup tig number : {}'.format(allOrigingroup_df.shape[0]))
        SetLog.error_out('grouptig2df key number : {}'.format(len(grouptig2df)))
    return allOrigingroup_df, grouptig2df

def read_file(movestat_xlsx):
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
    Wrong_group = []
    for num in range(1, 81):
        if num != 81:
            Oridf_query_df = Allorigroup[Allorigroup['Group'] == 'group{}'.format(num)]
            # print(Oridf_query_df.shape[0])
            # print(Oridf_query_df)
            # Oridf_query_df.to_excel('group73.origin_df.xlsx',header=False,index=False)
            # print('Ori'+'='*20)
            # print(num)
            Remove_query_df = movestat_df[movestat_df['pos1'] == 'group{}'.format(num)]#[0].unique().tolist()  # TODO 这里unique只保留靠前的group
            # print(Remove_query_df)
            # Remove_query_df.to_excel('group73.remove_df.xlsx',header=False,index=False)
            Remove_query_ctg = Remove_query_df['tigID'].tolist()
            # print(len(Remove_query_ctg))

            Add_query_df = movestat_df[movestat_df['pos2'] == 'group{}'.format(num)]
            # print(Add_query_df)
            # Add_query_df.to_excel('group73.add_df.xlsx',header=False,index=False)
            Add_query_ctg = Add_query_df['tigID'].tolist()
            # print(len(Add_query_ctg))

            # tmp = needmove_df[needmove_df[1] == 'group{}'.format(num)][0].tolist()
            # print(len(tmp))
            concat_list = []
            removed_query_df = Oridf_query_df[~Oridf_query_df['Contig'].isin(Remove_query_ctg)]
            # print(removed_query_df.shape[0])
            # removed_query_df.to_excel('group73.removed_df.xlsx',header=False,index=False)
            needadd_query_df = Allorigroup[Allorigroup['Contig'].isin(Add_query_ctg)]
            # print(needadd_query_df)
            result_query_df = pd.concat([removed_query_df,needadd_query_df],ignore_index=True)
            # result_query_df.to_excel('group73.result_df.xlsx',header=False,index=False)
            # print(removed_query_df)
            # print(result_query_df)
            result_query_df = result_query_df.iloc[:,[0,1,2]]
            result_query_df.rename(columns={'Contig':'#Contig'},inplace=True)
            # print(result_query_df)
            if Oridf_query_df.shape[0] - len(Remove_query_ctg) + len(Add_query_ctg) == result_query_df.shape[0]:
                SetLog.info_out('Group{} || The number of ctg moved is consistent with the number of ctg counted ~'.format(num))
                if output_path:
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
                    result_query_df.to_csv('{}/group{}.txt'.format(output_path,num),sep='\t',header=True,index=False)
        # break

    #     AddID_query = needmove_df[needmove_df[2] == 'group{}'.format(num)]
    #     # print(AddID_query)
    #     AddID_query = AddID_query[0].unique().tolist()
    #     Adddf_query = Allorigroup[Allorigroup['Contig'].isin(AddID_query)]
    #     Adddf_query = Adddf_query.drop_duplicates('Contig',keep='first')  #TODO 这里是对原始group里提取的adddf第一列去冗余，保留第一次出现的
    #     if len(AddID_query) == Adddf_query.shape[0]:
    #         SetLog.info_out(
    #             'The number of add rows extracted is consistent with the total number of add rows | [check ok ~]')
    #         # print(Adddf_query)
    #     else:
    #         SetLog.error_out('add contigID num != add df shape !!!')
    #         print('add contigID num : {}'.format(len(AddID_query)))
    #         print('add df rows : {}'.format(Adddf_query.shape[0]))
    #
    #     Resultdf_query = Oridf_query[~Oridf_query['Contig'].isin(RemoveID_query)]
    #     # print(Resultdf_query.shape[0])
    #     Resultdf_query = pd.concat([Resultdf_query,Adddf_query],ignore_index=True)
    #     # print(Resultdf_query)
    #     Resultdf_query.rename(columns={'Contig': '#Contig'}, inplace=True)
    #     Resultdf_query = Resultdf_query.iloc[:, 0:3]
    #
    #     if Oridf_query.shape[0] - len(RemoveID_query) + len(AddID_query) == Resultdf_query.shape[0]:
    #         SetLog.info_out('Group{} | The modified number of contigs is consistent with the calculated | [check ok ~~]'.format(num))
    #         Resultdf_query.to_csv('{}/group{}.txt'.format(output_path,num),sep='\t',header=True,index=False)
    #     else:
    #         SetLog.error_out('Group{} | The modified number of contigs is inconsistent with the calculated !!!'.format(num))
    #         SetLog.debug_out('Origin rows : {}'.format(Oridf_query.shape[0]))
    #         SetLog.debug_out('Remove contig num : {}'.format(len(RemoveID_query)))
    #         SetLog.debug_out('Add contig num : {}'.format(len(AddID_query)))
    #         SetLog.debug_out('Resultdf rows : {}'.format(Resultdf_query.shape[0]))
    #         Resultdf_query.to_csv('{}/group{}.txt'.format(output_path, num), sep='\t', header=True, index=False)
    #         Wrong_group.append(num)
    #     # print(num)
    #     break
    # print(Wrong_group)



        ###-------------Test duolicates-------------
        # tmp = [];dup_contig = []
        # for index,row in needmove_df[needmove_df[1] == 'group{}'.format(num)].iterrows():
        #     if row[0] not in tmp:
        #         tmp.append(row[0])
        #     else:
        #         dup_contig.append(row[0])
        # print(dup_contig)
        # break
        ##------------------------------------------

    # print(RemoveID_query)
    # print(needmove_df[needmove_df[1] == 'group{}'.format(num)])
    # print(len(RemoveID_query))
    # AddID_query = needmove_df[needmove_df[2] == 'group{}'.format(num)]
    # print(AddID_query)
    # # if AddID_query.shape[0] > 0: #TODO chr1做完后删除这一行
    # AddID_query = AddID_query[0].unique().tolist()
    # Adddf_query = Allorigroup[Allorigroup['Contig'].isin(AddID_query)]
    # if len(AddID_query) == Adddf_query.shape[0]:
    # SetLog.info_out('The number of add rows extracted is consistent with the total number of add rows | [check ok ~]')
    # print(Adddf_query)
    # else:
    # SetLog.error_out('add contigID num != add df shape !!!')
    # print('add contig num : {}'.format(len(AddID_query)))
    # print('add df rows : {}'.format(Adddf_query.shape[0]))
    # # Adddf_query.rename(columns={0:'Contig',1:'RECounts',2:'Length'},inplace=True)
    # # Adddf_query['Group'] = 'group{}'.format(num)

    # Removeddf_query = Oridf_query[~Oridf_query['Contig'].isin(RemoveID_query)]
    # if Oridf_query.shape[0] - len(RemoveID_query) == Removeddf_query.shape[0]:
    # SetLog.info_out('Origin_num - Remove_num == Removed_num | [check ok ~]')
    # print(Removeddf_query)
    # print('group{}:'.format(num))
    # print('Origin num : {}'.format(Oridf_query.shape[0]))
    # print('Remove num : {}'.format(len(RemoveID_query)))
    # print('Removed num : {}'.format(Removeddf_query.shape[0]))
    # Result_query = pd.concat([Removeddf_query,Adddf_query],ignore_index=True)
    # print(Result_query)
    # if Removeddf_query.shape[0] + Adddf_query.shape[0] == Result_query.shape[0]:
    # SetLog.info_out('Removed_num + Added_num  == Result_num | [check ok ~]')
    # print('Added num : {}'.format(Adddf_query.shape[0]))
    # print('Result num : {}'.format(Result_query.shape[0]))
    # Result_query.rename(columns={'Contig':'#Contig'},inplace=True)
    # Result_query = Result_query.iloc[:,0:3]
    # Result_query.to_csv('{}/group{}.txt'.format(output_path,num),sep='\t',index=False,header=True)
    # else:
    # SetLog.error_out('Removed_num + Added_num  != Result_num !')
    # print('Added num : {}'.format(Adddf_query.shape[0]))
    # print('Result num : {}'.format(Result_query.shape[0]))
    # else:
    # SetLog.error_out('Origin_num - Remove_num != result_num')
    # print('group{}:'.format(num))
    # print('Origin num : {}'.format(Oridf_query.shape[0]))
    # print('Remove num : {}'.format(len(RemoveID_query)))
    # print('Removed num : {}'.format(Removeddf_query.shape[0]))
    # print(needmove_df[needmove_df[1] == 'group{}'.format(num)])
    # print(needmove_df[needmove_df[1] == 'group{}'.format(num)][0].duplicated())


if __name__ == '__main__':
    # allorigingroup_df,grouptig2df = read_origin_group(r'E:\Badila\synteny_analysis\Badila染色体第二次校正\changeID_groups')
    allorigingroup_df,grouptig2df = read_origin_group(r'E:\Badila\synteny_analysis\Badila染色体第一次校正\矫正后\movetigstat\origin_group_test') #最初的 一次也没矫正过的
    # allorigingroup_df.to_excel(r'allorigingroup_df.xlsx',header=False,index=False)
    # allorigingroup_df,grouptig2df = read_origin_group(r'E:\Badila\synteny_analysis\Badila染色体第一次校正\矫正后\movetigstat\firstcorrected_group')
    # print(allorigingroup_df)
    # needmove_df = need2movetig(r'E:\Badila\synteny_analysis\Badila染色体第二次校正\矫正后\ZG_correct_round2.tigmovestat.allchr.xlsx',10)
    movestat_df = read_file(r'E:\Badila\synteny_analysis\Badila染色体第一次校正\矫正后\movetigstat\movestat.clean.xlsx')
    Reorganize_group(allorigingroup_df,movestat_df,r'E:\Badila\synteny_analysis\Badila染色体第一次校正\矫正后\movetigstat\firstcorrected_groupv0622')
    # Reorganize_group(allorigingroup_df,needmove_df,r'E:\Badila\synteny_analysis\Badila染色体第二次校正\矫正后\corrected_groups_round2')

    # allorigingroup_df, grouptig2df = read_origin_group(sys.argv[1])
    # print(allorigingroup_df)
    # needmove_df = need2movetig(sys.argv[2], sys.argv[3])
    # Reorganize_group(allorigingroup_df, needmove_df)
