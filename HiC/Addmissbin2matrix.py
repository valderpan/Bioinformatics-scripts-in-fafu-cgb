#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan

'''
Add the missing bin in the matirx file according to the output of Hic-Pro
'''

import pandas as pd
import sys
sys.path.append('/public1/home/stu_panhaoran/scripts')
import Fontcolor as Font

def add_missing_lines(bed_file,matrix_file,output_file):
    '''
    处理matirx矩阵文件中两个bin之间的互作值，添加未有的bin并赋值互作为0
    '''
    beddf = pd.read_table(bed_file, sep='\t', header=None)
    bin_numbers = beddf.iloc[beddf.shape[0] - 1,][3]

    df = pd.read_table(matrix_file,header=None,sep='\t',names=['firstbin','secondbin','frequency'],float_precision='round_trip')
    concat_list = []
    for i in df['firstbin'].unique().tolist():
        dd = df[df['firstbin'] == i]
        stars = i;ends=int(bin_numbers)
        append_list = []
        for j in range(stars,ends+1):
            if j not in dd['secondbin'].tolist():
                row = {'firstbin':i,'secondbin':j,'frequency':0}
                append_list.append(row)
        ds = dd.append(pd.DataFrame(append_list),ignore_index=True)
        ds = ds.sort_values(by='secondbin')
        concat_list.append(ds)
    result = pd.concat(concat_list,ignore_index=True)
    print(result)
    result.to_csv(output_file,sep='\t',index=False,header=True)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        Font.Scripts_tip('Add the missing bin in the matirx file according to the output of Hic-Pro')
        print('Usage:')
        print('\tpython {0} <bed.file> <matrix.file> <Correct_matrix.file>'.format(sys.argv[0]))
    else:
        add_missing_lines(sys.argv[1],sys.argv[2],sys.argv[3])
