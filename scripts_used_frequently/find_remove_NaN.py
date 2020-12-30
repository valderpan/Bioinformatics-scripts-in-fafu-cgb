#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2020/12/30

import sys
sys.path.append('/public1/home/stu_panhaoran/scripts')
import pandas as pd
import Fontcolor as Font

def parse_Nullrow(file):
    if file.endswith('.xlsx'):
        df = pd.read_excel(file,comment='#')
    elif file.endswith('.csv'):
        df = pd.read_csv(file,comment='#')
    else:
        df = pd.read_table(file,header=True,comment='#')

    nullrow = df[df.isnull().T.any()]
    nullrow.to_csv('Nullrows.tab',sep='\t',header=True,index=False)
    null_index = df[df.isnull().T.any()].index.tolist()
    df = df.drop(index=null_index,axis=1)
    df.to_csv('RemovedNull.tab',sep='\t',header=True,index=False)
    Font.Log_output('Number of rows with null values : {}'.format(len(null_index)))

if __name__ == '__main__':
    if len(sys.argv) != 2:
        Font.Scripts_tip('Find and remove the rows with Null values(Nan)')
        Font.Scripts_tip('\tNote:This script will create two files: <Nullrows.tab> <RemoveNull.tab>')
        print('Usage:')
        print('\tpython {0} <file>'.format(sys.argv[0]))
    else:
        parse_Nullrow(sys.argv[1])