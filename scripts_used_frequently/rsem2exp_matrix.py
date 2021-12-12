#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/12/12


import sys
sys.path.append('/share/home/stu_panhaoran/scripts')
import subprocess
from datetime import datetime
import pandas as pd
import argparse
import richlog as rl
from path import Path


__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20211212'
__version__ = 'v1.0'


def run_command(cmd):
    print(cmd)
    return_code = subprocess.call(cmd,shell=True)
    if return_code != 0 :
        print('ERROR:[{2}] Return code {0} when running the following command : {1} '.format(return_code, cmd,datetime.now()))


def get_RSEMdirs(dirpath):
    dirs = [dir for dir in Path(dirpath).dirs()]
    return dirs


def rm_dirsbam(dirs):
    for dir in dirs:
        run_command('rm {}/*.bam'.format(dir))


def gain_Samplematrix(dirs,Type):
    df_L = []
    for dir in dirs:
        ndir = dir.basename()
        df = pd.read_table('{}/RSEM.genes.results'.format(dir),sep='\t')
        if Type == 'FPKM':
            qdf = df.iloc[:,[0,-1]]
            qdf.rename(columns={'FPKM':ndir},inplace=True)
            df_L.append(qdf)
        elif Type == 'TPM':
            qdf = df.iloc[:,[0,-2]]
            qdf.rename(columns={'TPM': ndir}, inplace=True)
            df_L.append(qdf)
        else:
            rl.error('Error:Please enter the correct expression type: TPM or FPKM')
    return df_L


def merge2matrix(df_L,prefix,Type):
    mdf = df_L[0]
    for i in range(1,len(df_L)):
        mdf = pd.merge(mdf,df_L[i],on='gene_id')
    mdf.to_excel('{}.{}.xlsx'.format(prefix,Type),header=True,index=False)
    mdf.to_csv('{}.{}.csv'.format(prefix,Type),header=True,index=False)


def main(args):
    from rich.traceback import install
    install()
    dirspath = args.path
    Type = args.type
    oprefix = args.outprefix
    dirs = get_RSEMdirs(dirspath)
    rm_dirsbam(dirs)
    df_L = gain_Samplematrix(dirs,Type)
    merge2matrix(df_L,oprefix,Type)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''This script is used to convert RSEM quantitative results into a gene expression matrix''',
        usage="python {} -p rsem_out_path -t FPKM/TPM -p output_prefix ".format(sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,
                                                                            __mail__, __date__, __version__))
    parser.add_argument('-p', '--path', required=True, help='Input RSEM output path')
    parser.add_argument('-t', '--type', required=True,choices=['FPKM','TPM'], help='Input expression values type')
    parser.add_argument('-o', '--outprefix', required=True, help='Please specify the prefix of the output file')
    args = parser.parse_args()
    main(args)

