# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

# 计算表格中三列的平均值并将其显示在新增一列中
import pandas as pd
import os
import sys

os.chdir(r'D:\Result\test')
data = pd.read_excel('genes_fpkm.xls',sheet_name='spleen')
temp = data[['s35_FPKM','s46_FPKM','s57_FPKM']]
data['avgrage'] = temp.mean(axis=1)
data.to_excel('spleen_out.xlsx',index=False,header=True)
# print(data)