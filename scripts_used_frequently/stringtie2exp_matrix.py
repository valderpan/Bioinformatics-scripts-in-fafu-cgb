#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# @Time : 2020/10/11 11:04

'''
需求：Hisat2+stringtie表达量计算后将每个样本中基因的表达值整理为矩阵用于后续分析
'''

import pandas as pd
import re
import os
import sys

def get_fpkm_file(file_path):
    all_file_list = os.listdir(file_path)
    fpkm_file_list = [file for file in all_file_list if file.endswith('.tab')]
    return fpkm_file_list

def merge_file2matrix(file_path,file_list,output):
    df2merge_list = []
    for file in file_list:
        file_name = re.findall('([0-9A-Za-z\_\.\-]*)\.tab',file)[0]
        df = pd.read_table(file_path+'/'+file,sep='\t')
        df2merge = round(df.iloc[:,[0,-2]],2)
        df2merge.rename(columns={'FPKM':file_name},inplace=True)
        df2merge_list.append(df2merge)

    first_df2merge = df2merge_list[0]
    df2merge_dict = {}

    for i in range(1,len(df2merge_list)):
        df2merge_dict[i] = df2merge_list[i]
        if i < len(df2merge_list):
            first_df2merge = pd.merge(first_df2merge,df2merge_dict[i],on='Gene ID')

    first_df2merge.to_csv(output,header=True,index=False,sep='\t')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage:")
        print('\tpython {0} file_path output_file'.format(sys.argv[0]))
    else:
        file_list = get_fpkm_file(sys.argv[1])
        merge_file2matrix(sys.argv[1],file_list,sys.argv[2])