#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan

'''
需求：此脚本用来将每组等位中的三个/四个基因的cis-element转换为bed的格式用来画图，并统计每个gene上每个CRE出现的次数
'''

import pandas as pd
import re
import sys
import os
import argparse

__author__ = 'Haoran Pan'
__mail__ = 'haoran_pan@qq.com'
__date__ = '20201023'
__version__ = 'v1.0'


def gain_rowCRE_bed_allpick(input_path1,output_path1):
    '''
    统计每一行的3/4个等位基因内的所有CRE(包括与基因同向/反向的CRE)，生成bed文件
    '''
    inputpath = os.listdir(input_path1)
    filelist = [file for file in inputpath if file.endswith('.tab')]
    for file in filelist:
        df = pd.read_table(input_path1+'/'+file,sep='\t')
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
        rowCRE_filename = re.findall('([0-9A-ZAa-z\_]+)\.tab',file)[0]
        result.to_csv(output_path1+'/'+rowCRE_filename+'.bed',sep='\t',header=True,index=False)



def gain_rowCRE_bed_picksame(input_path2,output_path2):
    '''
    统计每一行的3/4个等位基因内的仅仅与基因同向的CRE，生成bed文件
    '''
    inputpath = os.listdir(input_path2)
    filelist = [file for file in inputpath if file.endswith('.tab')]
    for file in filelist:
        df = pd.read_table(input_path2+'/'+file,sep='\t')
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

            elif row['GeneID'][-1] == '+':
                if row['Strand'] == '+':
                    position_tmp = df.loc[[index], :]['Position'].tolist()[0][1:-1]
                    position_list = [int(i) for i in position_tmp.split(', ')]
                    for num in range(0, len(position_list)):
                        if num < len(position_list):
                            tmp = {'GeneID': row['GeneID'], 'CRE': row['Site Name'], 'Start': position_list[num],
                                   'End': position_list[num] + int(row['Length']), 'Strand': 'forward', 'direction': 1}
                            rowCRE_concat.append(pd.DataFrame(tmp, index=[0]))

        result = pd.concat(rowCRE_concat)
        rowCRE_filename = re.findall('([0-9A-ZAa-z\_\.]+)\.tab',file)[0]
        result.to_csv(output_path2+'/'+rowCRE_filename+'.bed',sep='\t',header=True,index=False)

def gain_rowCRE_info(input_path3,output_path3):
    '''
    将CRE与gene同向的基因进行分析，当四组等位在每一个CRE中均有分布时，且总的CRE分布统计数>总列数半数时输出统计情况
    '''

    inputpath = os.listdir(input_path3)
    filelist = [file for file in inputpath if file.endswith('.bed')]

    for file in filelist:
        df = pd.read_table(input_path3+'/'+file,sep='\t')
        df = df.groupby(['GeneID'])['CRE'].value_counts().unstack().reset_index()
        all_value = 0
        for i in range(1,df.shape[1]):
            # print([i for i in df.iloc[:,i].tolist() if i > 0 ])
            col_length = len([i for i in df.iloc[:,i].tolist() if i > 0 ])
            if col_length == df.shape[0]:
                all_value += 1
        if all_value >= df.shape[1] * 0.8:
            rowCRE_bedname = re.findall('([0-9a-zA-Z\_\.]+)\.bed',file)[0]
            df.to_csv(output_path3+'/'+rowCRE_bedname+'.info',sep='\t',header=True,index=False)


def main(args):
    choose = args.choose
    inputpath1 =args.inputpath1
    inputpath2 =args.inputpath2
    outputpath1 = args.output1
    outputpath2 = args.output2
    outputpath3 = args.output3
    if choose == 'all':
        gain_rowCRE_bed_allpick(inputpath1,outputpath1)
        gain_rowCRE_info(outputpath1, outputpath3)
    elif choose == 'same':
        gain_rowCRE_bed_picksame(inputpath1,outputpath2)
        gain_rowCRE_info(outputpath2, outputpath3)
    elif choose == 'none' :
        gain_rowCRE_info(inputpath2,outputpath3)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description=
        '''This script is used to convert the cis-element of the four genes in each group of alleles into bed format for drawing, and the output of the cis-element distribution on 3/4 alleles is similar (80%)
        If -c you choose "all" , use python {0} -c all -i1 -o1 -o3
        If -c you choose "same" , use python {1} -c same -i1 -o2 -o3
        If -c you choose "none" , use python {2} -c none -i2 -o3'''.format(sys.argv[0],sys.argv[0],sys.argv[0]),
        epilog= 'author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
    parser.add_argument('-c', '--choose',help='Select CRE in the same direction as the gene direction (-c you should choose "all"); '
                                              'Only select CRE in the same direction as the gene (-c you should choose "same"); '
                                              'just perform info on the output bed file Statistics (-c you should choose "none")')
    parser.add_argument('-i1', '--inputpath1',  help='Input the path where each row tab file is located')
    parser.add_argument('-i2', '--inputpath2',  help='Input the path where each row bed file is located')
    parser.add_argument('-o1', '--output1',  help='bed file output path(If -c you choose "all")')
    parser.add_argument('-o2', '--output2',  help='bed file output path(If -c you choose "same")')
    parser.add_argument('-o3', '--output3',  help='info file output path')
    args = parser.parse_args()
    main(args)