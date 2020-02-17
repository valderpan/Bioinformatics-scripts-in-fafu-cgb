# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

#遍历所提供的基因列表，在另一表格中将含有这一基因的一行输出并拼接为一个矩阵
import os
import pandas as pd

os.chdir(r'D:\Result\Transit\fpkm_caculate\单套')

df1 = pd.read_excel(r'C:\Users\dell\Desktop\PolyA\fpkm(单套).xlsx',sheet_name = 'ses-ID',header=None)
df2 = pd.read_excel('09.day_light_rhythm_2017.05.SES.fpkm_dataframe.xlsx')

# 生成含有目标geneID列所在位置索引的列表
number_list = []
for n in range(1,85):
    if n%2  == 1 :
        number_list.append(n)
# print(number_list)

#遍历有geneID的每一列
for num in number_list:
    #遍历每一列中的每一行，放到列表里
    gene_id_list=[]
    for i in df1[num].dropna():
        gene_id_list.append(i)
    print(gene_id_list)
    #读取列表里的第一个geneID
    first_data = df2[df2['gene_id'].str.contains(gene_id_list[0])]
    # print(first_data)
    file_dict = {}

    for j in gene_id_list[1:]:
        file_dict[j] = df2[df2['gene_id'].str.contains(j)]
        first_data = pd.concat([first_data,file_dict[j]])

    first_data.to_excel(str(num)+'.xlsx',index=False,header=True)
    #     print(first_data)