#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2020/11/26

'''
将目标基因的CRE提取出来转换为bed格式用于GSDS作图
'''
import sys
sys.path.append('../')
import pandas as pd
import Fontcolor as Font

def genebackground(geneIDfile):
    '''
    制作目标基因的背景文件
    :param geneIDfile:
    :return:
    '''
    df = pd.read_table(geneIDfile,names=['geneID'])
    df['Start'] = 0
    df['End'] = 2000
    df['Type'] = 'CDS'
    df.to_csv('genebackground.bed',sep='\t',header=None,index=False)


def gain_GeneBed_baseGeneID(geneIDfile,allgenomeCRE,output=0):
    '''
    将目标基因的CRE提取出来后转换为bed格式，即GeneID,CREstart,CREend,CREname
    :param geneIDfile:
    :param allgenomeCRE:
    :return:
    '''
    if allgenomeCRE.endswith('.xlsx'):
        allgeneCRE = pd.read_excel(allgenomeCRE)
    else:
        allgeneCRE = pd.read_table(allgenomeCRE,sep='\t')
    gene_df = pd.read_table(geneIDfile, header=None)
    gene_list = gene_df.iloc[:, 0].tolist()
    query_genedf = allgeneCRE[allgeneCRE['GeneID'].str[:-1].isin(gene_list)]
    query_genedf = query_genedf.iloc[:, [0,1,3,4,-1]]
    dd = query_genedf[query_genedf['GeneID'].str[-1] == query_genedf['Strand']]
    print(dd)
    rowCRE_concat = []
    for index,row in dd.iterrows():
        position_tmp = dd.loc[[index], :]['Position'].tolist()[0][1:-1]
        position_list = [int(i) for i in position_tmp.split(', ')]
        for num in range(0, len(position_list)):
            if num < len(position_list):
                tmp = {'GeneID': row['GeneID'], 'CRE': row['Site Name'], 'Start': position_list[num],
                       'End': position_list[num] + int(row['Length']), 'Strand': 'forward', 'direction': 1}
                rowCRE_concat.append(pd.DataFrame(tmp, index=[0]))
    result =pd.concat(rowCRE_concat)
    result['GeneID'] = result['GeneID'].str[:-1]
    result = result.iloc[:,[0,2,3,1]]
    print(result)
    if output:
        result.to_csv(output, sep='\t', header=None, index=False)
    return result


if __name__ == '__main__':
    if len(sys.argv) != 4:
        Font.Scripts_tip('The CRE of the query gene was extracted and converted to bed format for drawing with GSDS ')
        print('Usage:')
        print('\tpython {0} <querygeneIDlist> <allgenomeCRE> <Output.bed>'.format(sys.argv[0]))
    else:
        genebackground(sys.argv[1])
        gain_GeneBed_baseGeneID(sys.argv[1],sys.argv[2],sys.argv[3])