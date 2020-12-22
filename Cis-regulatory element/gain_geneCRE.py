#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2020/11/6

'''
找到目标基因的CRE并输出
'''


import pandas as pd
from plotnine import *
import sys
sys.path.append('../')
import Fontcolor as Font

# pd.set_option('display.max_columns',15)
def get_geneCRE_base_geneID(geneID_file,allgenomeCRE,output=0):
    allgeneCRE = pd.read_table(allgenomeCRE,sep='\t')
    querygeneID = pd.read_table(geneID_file,header=None).iloc[:,0].tolist()

    querygeneCRE = allgeneCRE[allgeneCRE['GeneID'].str[:-1].isin(querygeneID)]
    querygeneCRE = querygeneCRE[querygeneCRE['GeneID'].str[-1] == querygeneCRE['Strand']] #筛选同向的CRE
    # print(querygeneCRE)
    if output:
        # querygenefilename = re.findall('([0-9a-zA-Z\.]+)\.txt',geneID_file)[0]
        querygeneCRE.to_csv(output,sep='\t',header=True,index=False)
    return querygeneCRE

if __name__ == '__main__':
    if len(sys.argv) != 4:
        Font.Scripts_tip('Get the CRE of the query gene')
        print('Usage:')
        print('\tpython {0} <querygeneIDlist> <allgenomeCRE> <Output.tab>'.format(sys.argv[0]))
    else:
        get_geneCRE_base_geneID(sys.argv[1],sys.argv[2],sys.argv[3])