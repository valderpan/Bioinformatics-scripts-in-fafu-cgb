#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# @Time : 2021/08/09


import pandas as pd
import re
from path import Path
import sys
sys.path.append('/share/home/stu_panhaoran/scripts')
import Fontcolor as Font

def getFpkmfile(file_path):
    files = [file for file in Path(file_path).files() if file.endswith('.tab')]
    return files


def merge2matrix(files,output_file,valuetype):
    mergedf_list = []
    for file in files:
        file_name = re.findall('([0-9A-Za-z\_\.\-]*)\.tab', file)[0]
        df = pd.read_table(file,sep='\t',float_precision='round_trip')
        df2merge = df.loc[:,['Gene ID',valuetype]]
        df2merge['value'] = round(df2merge[valuetype],2)
        df2merge = df2merge.loc[:,['Gene ID','value']]
        df2merge.rename(columns={'value':file_name},inplace=True)
        mergedf_list.append(df2merge)

    df_1 = mergedf_list[0].sort_values(by='Gene ID')

    for i in range(1,len(mergedf_list)):
        df_1 = pd.merge(df_1,mergedf_list[i],on="Gene ID")

    df_1.to_csv(output_file,sep='\t',header=True,index=False)


if __name__ == '__main__':
    if len(sys.argv) == 4:
        files = getFpkmfile(sys.argv[1])
        merge2matrix(files,sys.argv[2],sys.argv[3])
    else:
        Font.Scripts_tip('Combine the expression values of each sample into a matrix after calculating the expressions from Hisat2+stringtie pipeline')
        print('Usage:')
        print('\tpython {0} [Path of tab files] <output_tab_file> <FPKM|TPM>'.format(sys.argv[0]))
