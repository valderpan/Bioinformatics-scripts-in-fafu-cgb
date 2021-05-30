#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/5/29

'''
%prog [tab files path] [output_path='./']

Use this script to convert all .tab files in the directory to .xlsx files
'''

import sys
import pandas as pd
from path import Path
import re

def get_tab_files(tab_path):
    files = [file for file in Path(tab_path).files() if file.endswith('.tab')]
    return files


def convert_format(files_list,output_path=0):
    if output_path:
        for file in files_list:
            file_prefix =  re.findall('([0-9a-zA-Z\.\_\-]+)\.tab',file.basename())[0]
            df = pd.read_table(file,sep='\t',names=['V1','V2','V3','V4','V5','V6','V7'])
            df.to_excel('{}/{}.xlsx'.format(output_path,file_prefix),header=False,index=False)
    else:
        for file in files_list:
            file_prefix =  re.findall('([0-9a-zA-Z\.\_\-]+)\.tab',file.basename())[0]
            print(file_prefix)
            df = pd.read_table(file,sep='\t',names = ['V1','V2','V3','V4','V5','V6','V7'])
            df.to_excel('{}.xlsx'.format(file_prefix),header=False,index=False)


if __name__ == '__main__':
    from optparse import OptionParser

    p = OptionParser(__doc__)
    opts, args = p.parse_args()
    if len(args) == 2:
        tab_files_path = args[0]
        output_file_path = args[1]
        tabfiles_list = get_tab_files(tab_files_path)
        convert_format(tabfiles_list,output_file_path)
    else:
        sys.exit(p.print_help())
