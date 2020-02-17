# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

# 根据某一列数值计算log2（将此列数值+1后计算log2），并在新的一列输出
import pandas as pd
import os
import math
import sys

os.chdir(r"C:\Users\dell\Desktop\Goldfish\data")
df = pd.read_csv('testis_fpkm.csv',header = None)
df.columns = ['Gene_Id','gene_id_fpkm','Homologs in goldfish','Homo_fpkm']
log2FC_gene = open('log2FC_gene.xlsx','w')

lon2FC_homo = open('log2FC_homo.xlsx','w')
temp = df[['gene_id_fpkm','Homo_fpkm']]

for i in df['gene_id_fpkm']:
    i +=1
    i = math.log2(i)
    i = str(i)
    log2FC_gene.write(i+'\n')
log2FC_gene.close()

for j in df['Homo_fpkm']:
    j +=1
    j = math.log2(j)
    j = str(j)
    lon2FC_homo.write(j+'\n')
log2FC_gene.close()