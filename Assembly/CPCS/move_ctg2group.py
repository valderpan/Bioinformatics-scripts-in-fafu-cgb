#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/3/18

'''
根据ROC22_6的输出结果进行ctg的转移
'''

import pandas as pd
import CPCS3
from path import Path
import shutil

def moveCtg2group(nonDict,chrinfo,oldpath,outpath):
    # for key in nonDict.keys():
    with open(chrinfo) as f:
        Connum = []
        Nonnum = []
        lines = (line.strip() for line in f)
        for line in lines:
            if line == 'end':
                continue
            line_list = line.split(' need move to ')
            # print(line_list)
            if line_list[0].split('-')[0] not in Nonnum and line_list[1].split('-')[0] not in Connum:
                Nonnum.append(line_list[0].split('-')[0])
                Connum.append(line_list[1].split('-')[0])
                nondf = pd.read_table(oldpath + '\\' + 'group' + line_list[0].split('-')[0] + '.txt', sep='\t')
                condf = pd.read_table(oldpath + '\\' + 'group' + line_list[1].split('-')[0] + '.txt', sep='\t')
                moveCtg = nonDict[line_list[0]]['ROC_chr'].unique().tolist()
                #print(moveCtg)

                # print(condf)

                # moveCtg_df = nonDict[line_list[0]][nonDict[line_list[0]]['ROC_chr'].isin(moveCtg)]
                moveCtg_df = nondf[nondf['#Contig'].isin(moveCtg)]
                condf = pd.concat([condf,moveCtg_df])
                condf.to_csv(outpath+'\\'+'group'+line_list[1].split('-')[0]+'.txt',sep='\t',header=True,index=False)
                # print(condf)
                # print(condf.shape[0])
                for i in moveCtg:
                    nondf = nondf[~nondf['#Contig'].isin([i])]
                # print(df)
                nondf.to_csv(outpath+'\\'+'group'+line_list[0].split('-')[0]+'.txt',sep='\t',header=True,index=False)
            elif line_list[0].split('-')[0] in Nonnum and line_list[1].split('-')[0] not in Connum:
                Connum.append(line_list[1].split('-')[0])
                nondf = pd.read_table(outpath + '\\' + 'group' + line_list[0].split('-')[0] + '.txt', sep='\t')
                condf = pd.read_table(oldpath + '\\' + 'group' + line_list[1].split('-')[0] + '.txt', sep='\t')

                moveCtg = nonDict[line_list[0]]['ROC_chr'].unique().tolist()
                moveCtg_df = nondf[nondf['#Contig'].isin(moveCtg)]
                condf = pd.concat([condf, moveCtg_df])
                condf.to_csv(outpath + '\\' + 'group' + line_list[1].split('-')[0] + '.txt', sep='\t', header=True,index=False)
                # print(condf)
                # print(condf.shape[0])
                for i in moveCtg:
                    nondf = nondf[~nondf['#Contig'].isin([i])]
                # print(df)
                nondf.to_csv(outpath + '\\' + 'group' + line_list[0].split('-')[0] + '.txt', sep='\t', header=True,index=False)

            elif line_list[0].split('-')[0] not in Nonnum and line_list[1].split('-')[0] in Connum:
                Nonnum.append(line_list[0].split('-')[0])
                nondf = pd.read_table(oldpath + '\\' + 'group' + line_list[0].split('-')[0] + '.txt', sep='\t')
                condf = pd.read_table(outpath + '\\' + 'group' + line_list[1].split('-')[0] + '.txt', sep='\t')
                moveCtg = nonDict[line_list[0]]['ROC_chr'].unique().tolist()
                #print(moveCtg)

                # print(condf)

                # moveCtg_df = nonDict[line_list[0]][nonDict[line_list[0]]['ROC_chr'].isin(moveCtg)]
                moveCtg_df = nondf[nondf['#Contig'].isin(moveCtg)]
                condf = pd.concat([condf,moveCtg_df])
                condf.to_csv(outpath+'\\'+'group'+line_list[1].split('-')[0]+'.txt',sep='\t',header=True,index=False)
                # print(condf)
                # print(condf.shape[0])
                for i in moveCtg:
                    nondf = nondf[~nondf['#Contig'].isin([i])]
                # print(df)
                nondf.to_csv(outpath+'\\'+'group'+line_list[0].split('-')[0]+'.txt',sep='\t',header=True,index=False)
            elif line_list[0].split('-')[0] in Nonnum and line_list[1].split('-')[0] in Connum:
                nondf = pd.read_table(outpath + '\\' + 'group' + line_list[0].split('-')[0] + '.txt', sep='\t')
                condf = pd.read_table(outpath + '\\' + 'group' + line_list[1].split('-')[0] + '.txt', sep='\t')
                moveCtg = nonDict[line_list[0]]['ROC_chr'].unique().tolist()
                # print(moveCtg)

                # print(condf)

                # moveCtg_df = nonDict[line_list[0]][nonDict[line_list[0]]['ROC_chr'].isin(moveCtg)]
                moveCtg_df = nondf[nondf['#Contig'].isin(moveCtg)]
                condf = pd.concat([condf, moveCtg_df])
                condf.to_csv(outpath + '\\' + 'group' + line_list[1].split('-')[0] + '.txt', sep='\t', header=True,index=False)
                # print(condf)
                # print(condf.shape[0])
                for i in moveCtg:
                    nondf = nondf[~nondf['#Contig'].isin([i])]
                # print(df)
                nondf.to_csv(outpath + '\\' + 'group' + line_list[0].split('-')[0] + '.txt', sep='\t', header=True,index=False)

def copy_nomovefile(oldpath,newpath):
    oldfiles = Path(oldpath).files()
    newfiles = Path(newpath).files()
    olds = [file.basename() for file in oldfiles]
    news = [file.basename() for file in newfiles]
    for old in olds:
        if old not in news:
            shutil.copy2(oldpath+'\\'+old,newpath)



if __name__ == '__main__':
    nondict = CPCS3.Store_all_nonblocks(r'E:\Badila\synteny_analysis\C03\parsed1')
    moveCtg2group(nondict,r'E:\Badila\synteny_analysis\C03\parsed1\C03.parse1.txt',r'E:\Badila\synteny_analysis\C03\new',r'E:\Badila\synteny_analysis\C03\parsed1\group_parsed1')
    copy_nomovefile(r'E:\Badila\synteny_analysis\C03\new',r'E:\Badila\synteny_analysis\C03\parsed1\group_parsed1')