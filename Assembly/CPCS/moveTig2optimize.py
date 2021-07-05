#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/6/11

'''
%prog <Origin_groupfile_path>

此脚本用来解完麻子后，根据movestat表将groupXX.txt中的congtig在各个group之间进行移动，并最后生成新的groupXX.txt
'''

import re
import sys
import argparse
import pandas as pd
from path import Path
from loguru import logger


__author__ = 'PHR'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20210705'
__version__ = 'v2.1'


logger.remove()
logger.add(sys.stdout,
           format="<green>{time:[YYYY-MM-DD HH:mm:ss]}</green>  | <level>{level}</level>  | <level>{message}</level>",
           colorize=True)


def info_out(input_str):
    logger.info(input_str)


def warning_out(input_str):
    logger.warning(input_str)


def error_out(input_str):
    logger.error(input_str)


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
        info_out('There is no duplicate tig in allOrigingroup dataframe! | [check ok ~]')
    else:
        error_out('AllOrigingroup dataframe contains duplicate tig!')
        error_out('Allorigingroup tig number : {}'.format(allOrigingroup_df.shape[0]))
        error_out('grouptig2df key number : {}'.format(len(grouptig2df)))
    return allOrigingroup_df, grouptig2df


def read_file(movestat_file):
    info_out('Read movestat file')
    if movestat_file.endswith('.xlsx'):
        movestat_df = pd.read_excel(movestat_file,names=['tigID','pos1','pos2'])
    else:
        movestat_df = pd.read_table(movestat_file, names=['tigID', 'pos1', 'pos2'],sep='\t')
    return movestat_df


def Reorganize_group(Allorigroup, movestat_df, output_path=0):
    '''
    对原有的group进行重排，重排后的groupxx.txt = origin - remove + add
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

            removed_query_df = Oridf_query_df[~Oridf_query_df['Contig'].isin(Remove_query_ctg)]

            needadd_query_df = Allorigroup[Allorigroup['Contig'].isin(Add_query_ctg)]

            result_query_df = pd.concat([removed_query_df,needadd_query_df],ignore_index=True)

            result_query_df = result_query_df.iloc[:,[0,1,2]]
            result_query_df.rename(columns={'Contig':'#Contig'},inplace=True)

            if Oridf_query_df.shape[0] - len(Remove_query_ctg) + len(Add_query_ctg) == result_query_df.shape[0]:
                info_out('Group{} || The number of ctg moved is consistent with the number of ctg counted ~'.format(num))
                if output_path:
                    result_query_df.drop_duplicates('#Contig',keep='first',inplace = True)
                    result_query_df.to_csv('{}/group{}.txt'.format(output_path,num),sep='\t',header=True,index=False)

            else:
                error_out('Group{} || The number of ctg moved is inconsistent with the number of ctg counted !!!'.format(num))
                error_out('Group{} Origin ctg : {}'.format(num,Oridf_query_df.shape[0]))
                error_out('Group{} need remove ctg : {}'.format(num,len(Remove_query_ctg)))
                error_out('Group{} need add ctg : {}'.format(num,len(Add_query_ctg)))
                error_out('Group{} result ctg should is {} but it is {}'.format(num,Oridf_query_df.shape[0] - len(Remove_query_ctg) + len(Add_query_ctg),result_query_df.shape[0]))
                for ctg in Remove_query_df['tigID'].tolist():
                    if ctg not in Oridf_query_df['Contig'].tolist():
                        warning_out('Group{} {} in Remove but it not in Origin'.format(num,ctg))
                if output_path:
                    result_query_df.drop_duplicates('#Contig', keep='first',inplace = True)
                    result_query_df.to_csv('{}/group{}.txt'.format(output_path,num),sep='\t',header=True,index=False)

def main(args):
    origin_path = args.originalpath
    input_file = args.movestat
    moved_path = args.movedpath
    allorigingroup_df, grouptig2df = read_origin_group(origin_path)
    movestat_df = read_file(input_file)
    Reorganize_group(allorigingroup_df,movestat_df,moved_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''This script is used to move the congtig in groupXX.txt between groups according to the movestat file after chromosome correction, and finally generate a new groupXX.txt''',
        usage="python {} -p original_path -i movestat_file -o output_path".format(sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__, __mail__, __date__,
                                                                            __version__))
    parser.add_argument('-p', '--originalpath', required=True, help='Inpout the path to the original groupXX.txt file')
    parser.add_argument('-i', '--movestat', required=True, help="Input contig's move stat file")
    parser.add_argument('-o', '--movedpath', required=True, help='Specify the path to the output file')
    args = parser.parse_args()
    main(args)

