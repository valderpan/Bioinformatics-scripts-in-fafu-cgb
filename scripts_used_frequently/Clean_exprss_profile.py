#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2020/11/27

'''
将合并的带有3个生物学重复的样本取均值构建表达谱
注意：运行前要根据生物学重复排好顺序,eg:L1-1,L1-2,L1-3,L2-1,L2-2,L2-3...
'''
import pandas as pd
import numpy as np
import re
import sys

def remove_duplicate_value_inlist(query_list):
    '''
    移除列表中重复的值
    '''
    new_list = []
    for i in query_list:
        if i not in new_list:
            new_list.append(i)
    return new_list


def average_rnaseq_expression(raw_expression_file,output_file=0):
    '''
    将表达谱中3个生物学重复取平均值并输出
    '''
    df = pd.read_excel(raw_expression_file,index_col=0)
    old_columns_names = [column for column in df]
    old_columns = [re.findall('([a-zA-Z0-9]+)\-',column)[0] for column in old_columns_names]  #TODO 这里要根据每个sample的ID进行手动修改
    old_columns = remove_duplicate_value_inlist(old_columns)
    df = df.groupby(np.arange(len(df.columns))//3,axis=1).mean()
    df = round(df, 2)
    df = df.reset_index()
    if output_file:
        df.to_excel(output_file,header=True,index=False)
    return df


def fpkm_filter(exp_profile,clean_exp_profile=0):
    '''
    将各个组织中的fpkm进行过滤，只要有两个组织中的fpkm>1,则保留这个gene
    :return:
    '''
    df = exp_profile
    concat_row = []
    for index,row in df.iterrows():
        num = 0
        for i in row.tolist()[1:]:
            if int(i)>1:
                num +=1
        if num >2 :
            # print(row)
            # print(df.loc[[index],:])
            concat_row.append(df.loc[[index],:])
    result = pd.concat(concat_row)
    if clean_exp_profile:
        result.to_excel(clean_exp_profile,header=True,index=False)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage:")
        print("\tpython {0} RawExpressionFile outputfile".format(sys.argv[0]))
    else:
        df = average_rnaseq_expression(sys.argv[1])
        fpkm_filter(df,sys.argv[2])
#
# df = average_rnaseq_expression(r'E:\\潘浩然\\Result\\Transit\\fpkm_caculate\\全套\\01.SES-208_2016.09-SingleEnd.fpkm_dataframe.xlsx')
# fpkm_filter(df,r'E:\潘浩然\Result\Transit\fpkm_caculate\cleaned\全套\SES_leaf.clean.xlsx')