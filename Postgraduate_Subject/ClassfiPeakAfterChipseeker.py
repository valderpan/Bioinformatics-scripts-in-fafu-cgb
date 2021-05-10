#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/4/29

'''
此脚本用于将chipseeker注释后的peaks再次进行分类，最后分成七类：Intron、Exon、Proximal Upstream、Intergenic、Proximal Downstream、3' UTR、5' UTR
'''
import pandas as pd

def classfipeakanno(peakanno,output):
    if peakanno.endswith('.xlsx'):
        df = pd.read_excel(peakanno)
    elif peakanno.endswith('.csv'):
        df = pd.read_csv(peakanno)
    else:
        df = pd.read_table(peakanno,sep='\t')
    new_anno = []
    for index,row in df.iterrows():
        if row['annotation'].startswith('Intron'):
            new_anno.append('Intron')
        elif row['annotation'].startswith('Exon'):
            new_anno.append('Exon')
        elif row['annotation'] == 'Distal Intergenic':
            new_anno.append('Intergenic')
        elif row['annotation'] == 'Promoter (<=1kb)':
            new_anno.append('Proximal Upstream')
        elif row['annotation'] == 'Promoter (1-2kb)':
            new_anno.append('Proximal Upstream')
        elif row['annotation'] == 'Promoter (2-3kb)':
            new_anno.append('Intergenic')
        elif row['annotation'] == 'Downstream (<1kb)':
            new_anno.append('Proximal Downstream')
        elif row['annotation'] == 'Downstream (1-2kb)':
            new_anno.append('Intergenic')
        elif row['annotation'] == 'Downstream (2-3kb)':
            new_anno.append('Intergenic')
        elif row['annotation'] == "3' UTR":
            new_anno.append("3' UTR")
        elif row['annotation'] == "5' UTR":
            new_anno.append("5' UTR")


    df['newanno'] = new_anno
    print(df)
    df.to_excel(output,header=True,index=False)

if __name__ == '__main__':
    classfipeakanno(r'E:\潘浩然\调控元件生信任务\analysis\Sugarcane\OsZmSb\Os_peak.txt',
                    r'E:\潘浩然\调控元件生信任务\analysis\Sugarcane\OsZmSb\Os_peak_mutated.xlsx')