# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

# 将Trinity得到每个表达量文件合并为一个表达矩阵

import os
import pandas as pd
import re

pd.set_option('display.max_columns', 10)
os.chdir(r'D:\Result\Transit\fpkm_caculate\test')

fpkm_file_list = os.listdir('./')
fpkm_file = []
for i in fpkm_file_list:
    j = re.findall('([a-zA-Z0-9\-]+)\.',i)
    fpkm_file.append(j[0])
# print(fpkm_file)

first_data = pd.read_table(fpkm_file_list[0],names=['gene_id',fpkm_file[0]])
first_data = first_data.drop(first_data.index[0])
file_dict = {}
# print(first_data)

for i in range(1,len(fpkm_file_list)):
    file_dict[i] = pd.read_table(fpkm_file_list[i],names=['gene_id',fpkm_file[i]])
    if i < len(fpkm_file_list):
        first_data = pd.merge(first_data,file_dict[i],on='gene_id')

first_data.to_excel('fpkm_dataframe.xlsx',header=True,index=False)
