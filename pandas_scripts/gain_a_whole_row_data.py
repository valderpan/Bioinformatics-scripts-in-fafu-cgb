# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

# 从一个表达矩阵中根据特定字符提取一张行数据
import os
import pandas as pd

os.chdir(r'D:\Result\Transit\fpkm_caculate')

data = pd.read_excel('09.day_light_rhythm_2017.05.SES.fpkm_dataframe.xlsx')
# print(data)
output = open('putput.xlsx','w')
# seq_list = []
with open('test.txt') as f:
    for line in f:
        line = line.rstrip()
        print(data[data.gene_id.str.contains(line)])
#         seq_list.append(line)
# # print(seq_list)
# output_list = []
# for i in seq_list:
#     gene_fpkm = data[data.gene_id.str.contains(i)]
#     output_list.append(gene_fpkm)
# print(output_list)
