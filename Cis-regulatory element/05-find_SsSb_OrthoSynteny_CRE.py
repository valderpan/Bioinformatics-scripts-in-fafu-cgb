#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# @Time : 2020/8/5 22:08

import pandas as pd
import sys
import os

with open(sys.argv[1]) as f:
    gene_ID = [line.strip() for line in f ]
df = pd.read_table(sys.argv[2],sep='\t')
print(df)
concat_list= []
for index,row in df.iterrows():
    if row['GeneID'][:-1] in gene_ID:
        #print(row['GeneID'][:-1])
        concat_list.append(df.loc[[index],:])

result = pd.concat(concat_list)
result.to_csv('SsSb_OrthoSynteny_CRE.tab',sep='\t',header=True,index=False)
