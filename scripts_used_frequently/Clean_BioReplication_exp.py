#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/09/05

'''
%prog <Raw_exp.xlsx> <sample nums> [Clean_exp.xlsx]

将合并的带有3/2个生物学重复的样本取均值构建表达谱
注意：运行前要根据生物学重复排好顺序且所有样本列数和为3/2的倍数,
eg:L_1-1,L_1-2,L_1-3,L_2-1,L_2-2,L_2-3...
'''

import re
import sys
import pandas as pd
import numpy as np


def remove_duplicate_value_inlist(query_list):
    '''
    移除列表中重复的值
    '''
    new_list = []
    for i in query_list:
        if i not in new_list:
            new_list.append(i)
    return new_list


def average_rnaseq_expression(raw_expression_file,sample_num,output_file):
    '''
    将表达谱中3个生物学重复取平均值并输出
    '''
    df = pd.read_excel(raw_expression_file,index_col=0)
    old_columns_names = [column for column in df]
    old_columns = [re.findall('([a-zA-Z0-9]+)\-',column)[0] for column in old_columns_names]
    old_columns = remove_duplicate_value_inlist(old_columns)
    df = df.groupby(np.arange(len(df.columns))//int(sample_num),axis=1).mean()
    df = round(df, 2)
    cols = [num for num in  range(0,df.shape[1])]
    newcols = dict(zip(cols,old_columns))
    df = df.reset_index()
    df.rename(columns=newcols,inplace=True)
    df.to_excel(output_file,header=True,index=False)



if __name__ == '__main__':
    from rich.traceback import install
    from optparse import OptionParser

    p = OptionParser(__doc__)
    opts, args = p.parse_args()
    if len(args) == 3:
        install()
        raw_file = args[0]
        samples = args[1]
        clean_file = args[2]
        average_rnaseq_expression(raw_file,samples,clean_file)
    else:
        sys.exit(p.print_help())
