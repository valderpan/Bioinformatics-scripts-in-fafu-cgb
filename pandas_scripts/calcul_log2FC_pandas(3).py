# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

#利用pandas对符合条件的行进行筛选
import pandas as pd
import os
import sys

os.chdir(r'C:\Users\dell\Desktop\Goldfish\data')
pd.set_option('display.max_columns', 10)
df = pd.read_excel('fpkm_log2FC.xlsx',sheet_name=sys.argv[1])

# 筛选的条件为两列log2值中有一列要>10，且最后相减得到的log2FC要>1或<-1
data = df[(df.gene_id_fpkm >=10) | (df.Homo_fpkm >=10)]
select_data = data[(data.log2FC >= 1) | (data.log2FC <= -1)]
# print(select_data)
select_data.to_excel(sys.argv[2],index=False,header=True)