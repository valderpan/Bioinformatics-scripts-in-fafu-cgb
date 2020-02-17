# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

# 讲xlsx文件储存在Dataframe中，并将两列中含有缺失值及congtig的行删除！
import os
import pandas as pd
os.chdir(r'C:\Users\dell\Desktop\Goldfish')
data = pd.read_excel('Supplementary Data 11 - Homolog gene pairs list for Subgenome A and B.xlsx')
# print(data[['Gene Id','Homologs in goldfish']])
data = data.astype(str) #必加！
data = data[~ data['Gene Id'].str.contains('nan') & ~ data['Gene Id'].str.contains('contig')]
data = data[~ data['Homologs in goldfish'].str.contains('nan') & ~ data['Homologs in goldfish'].str.contains('contig')]

data.to_excel('output.xlsx',index=False,header=True)