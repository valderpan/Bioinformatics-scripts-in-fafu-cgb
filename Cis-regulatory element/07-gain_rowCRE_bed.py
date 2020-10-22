#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan

'''
需求：此脚本用来将每组等位中的四个基因的cis-element转换为bed的格式用来画图
'''

import pandas as pd
import re
import sys

def gain_rowCRE_bed(rowCRE_file,output_path):
    df = pd.read_table(rowCRE_file,sep='\t')
    df = df.iloc[:,[0,1,2,3,4,7]]

    rowCRE_concat = []
    for index,row in df.iterrows():
        if row['GeneID'][-1] == '-':
            if row['Strand'] == '-':
                position_tmp = df.loc[[index],:]['Position'].tolist()[0][1:-1]
                position_list = [int(i) for i in position_tmp.split(', ')]
                for num in range(0,len(position_list)):
                    if num < len(position_list):
                        tmp = {'GeneID':row['GeneID'],'CRE':row['Site Name'],'Start':position_list[num],
                               'End':position_list[num]+int(row['Length']),'Strand':'forward','direction':1}
                        rowCRE_concat.append(pd.DataFrame(tmp, index=[0]))
            elif row['Strand'] == '+':
                position_tmp = df.loc[[index], :]['Position'].tolist()[0][1:-1]
                position_list = [int(i) for i in position_tmp.split(', ')]
                for num in range(0, len(position_list)):
                    if num < len(position_list):
                        tmp = {'GeneID': row['GeneID'], 'CRE': row['Site Name'], 'Start': position_list[num] - int(row['Length']),
                               'End': position_list[num], 'Strand': 'reverse', 'direction': -1}
                        rowCRE_concat.append(pd.DataFrame(tmp, index=[0]))
        elif row['GeneID'][-1] == '+':
            if row['Strand'] == '+':
                position_tmp = df.loc[[index], :]['Position'].tolist()[0][1:-1]
                position_list = [int(i) for i in position_tmp.split(', ')]
                for num in range(0, len(position_list)):
                    if num < len(position_list):
                        tmp = {'GeneID': row['GeneID'], 'CRE': row['Site Name'], 'Start': position_list[num],
                               'End': position_list[num] + int(row['Length']), 'Strand': 'forward', 'direction': 1}
                        rowCRE_concat.append(pd.DataFrame(tmp, index=[0]))
            elif row['Strand'] == '-':
                position_tmp = df.loc[[index], :]['Position'].tolist()[0][1:-1]
                position_list = [int(i) for i in position_tmp.split(', ')]
                for num in range(0, len(position_list)):
                    if num < len(position_list):
                        tmp = {'GeneID': row['GeneID'], 'CRE': row['Site Name'],
                               'Start': position_list[num] - int(row['Length']),
                               'End': position_list[num], 'Strand': 'reverse', 'direction': -1}
                        rowCRE_concat.append(pd.DataFrame(tmp, index=[0]))

    result = pd.concat(rowCRE_concat)
    rowCRE_filename = re.findall('([0-9A-ZAa-z\_]+)\.tab',rowCRE_file)[0]
    result.to_csv(output_path+'/'+rowCRE_filename+'.bed',sep='\t',header=True,index=False)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage:')
        print('\tpython {0} row.tab output_path'.format(sys.argv[0]))
    else:
        gain_rowCRE_bed(sys.argv[1],sys.argv[2])
