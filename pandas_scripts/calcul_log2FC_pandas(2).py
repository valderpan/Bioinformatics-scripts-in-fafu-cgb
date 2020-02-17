# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

# 将两列log2值相减并将结果输出在新的一列
import pandas as pd
import os
import sys

os.chdir(r'C:\Users\dell\Desktop\Goldfish\data')
pd.set_option('display.max_columns', 10)
df = pd.read_csv(sys.argv[1])
df.columns = ['Gene_Id','gene_id_fpkm','gene_log2FC','Homologs in goldfish','Homo_fpkm','Homo_log2FC']
df['log2FC'] = (df['Homo_log2FC']) - (df['gene_log2FC'])
df.to_excel(sys.argv[2],index=False,header=True)

