#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/9/9

'''
%prog <KaKsfile_path> <output.tab>

Extract the Ks values in each sample and concat
'''

import re
import sys
from path import Path
import pandas as pd

def Get_files(file_path):
    files = [i for i in Path(file_path).files() if i.endswith('.KaKs.result')]
    return files

def Read_files(files):
    C_L = []
    for file in files:
        name = re.findall('([0-9a-zA-Z]+)\.KaKs',file)[0]
        df = pd.read_table(file,sep='\t')
        # print(df)
        df['name'] = name
        df = df[df['dS-ng'] > 0]
        qdf = df.iloc[:,[0,3]]
        C_L.append(qdf)
    return  C_L


def ConcatOutput(c_l,output):
    result = pd.concat(c_l)
    result.to_csv(output,sep='\t',header=True,index=False)


if __name__ == "__main__":
    from optparse import OptionParser
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args)==2:
        file_path = args[0]
        output = args[1]
        files = Get_files(file_path)
        c_l = Read_files(files)
        ConcatOutput(c_l,output)
    else:
        sys.exit(p.print_help())
