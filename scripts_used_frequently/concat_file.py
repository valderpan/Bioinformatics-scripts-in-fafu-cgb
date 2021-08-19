#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/8/18

'''
%prog <file_path> <suffix> <output.xlsx>

Concat files under the path
'''

import sys
import pandas as pd
from path import Path


def get_files(file_path,suffix):
    files = [i for i in Path(file_path).files() if i.endswith(suffix)]
    return files


def concat_files(files):
    L = []
    for file in files:
        if file.endswith('.csv'):
            df = pd.read_csv(file)
            L.append(df)
        elif file.endswith('.xlsx'):
            df = pd.read_excel(file)
            L.append(df)
        else:
            df = pd.read_table(file)
            L.append(df)
    result = pd.concat(L)
    return result


def output(result,output):
    result.to_excel(output,header=True,index=False)


if __name__ == "__main__":
    from optparse import OptionParser
    from rich.traceback import install
    install()
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args)==3:
        file_path = args[0]
        suffix = args[1]
        Output = args[2]
        files = get_files(file_path,suffix)
        result = concat_files(files)
        output(result,Output)
    else:
        sys.exit(p.print_help())

