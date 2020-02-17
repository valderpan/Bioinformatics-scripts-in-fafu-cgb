# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

'''
此脚本用于根据DataframeA中的一列数，
从DataframeB中找到相应的列，并把其对应的值打印出来
'''
import pandas as pd
import os


pd.set_option('display.max_columns', 4)
os.chdir(r'C:\Users\dell\Desktop\Goldfish')
raw_data= pd.read_excel('Supplementary Data 11 - Homolog gene pairs list for Subgenome A and B.xlsx')
# 将header中的空格替换为下划线，方便于后面操作
raw_data.columns = [c.replace(' ','_') for c in raw_data.columns]

# 删除这两列中包含有NA和contig的行
Filter_data = raw_data.dropna(subset = ['Gene_Id'])
Filter_data = Filter_data.dropna(subset=['Homologs_in_goldfish'])
Filter_data = Filter_data[~Filter_data['Gene_Id'].str.contains('contig')]
Filter_data = Filter_data[~Filter_data['Homologs_in_goldfish'].str.contains('contig')]

# 截取列中字符串的长度
Filter_data['Gene_Id'] = Filter_data['Gene_Id'].str[0:14]
Filter_data['Homologs_in_goldfish'] = Filter_data['Homologs_in_goldfish'].str[0:14]

# print(Filter_data)

fpkm_data = pd.read_excel('genes_fpkm_1.xlsx',sheet_name= 'spleen')
# 将一行中的多个数值拆分为多行数值，此操作还不成熟，会使得列位置变化
fpkm_data = fpkm_data.drop('gene_short_name',axis=1).join(fpkm_data['gene_short_name'].str.split(',',expand=True).stack().reset_index(level=1,drop = True).rename('ID'))
# print(fpkm_data)
# 重新规定列的位置
fpkm_order = ['ID','avgrage']
fpkm_data = fpkm_data[fpkm_order]
# print(fpkm_data)

# 将Dataframe的每一行存为字典，ID列为key，average列为value
fpkm_dict = fpkm_data.set_index('ID').to_dict()['avgrage']

# 找到对应的值并打印，以下操作还不成熟，需要分布操作
for i in Filter_data['Gene_Id']:
    if i in fpkm_dict.keys():
        print(i,fpkm_dict[i])
    else:
        print(i)


for j in Filter_data['Homologs_in_goldfish']:
    if j in fpkm_dict.keys():
        print(j,fpkm_dict[j])
    else:
        print(j)