# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

import os
import pandas as pd

pd.set_option('display.max_columns', 15)
os.chdir(r'C:\Users\dell\Desktop\Goldfish\enrich')
### 以下是GO步骤
# df =pd.read_excel('Supplementary Data 12.xlsx',sheet_name='GO_list')
df =pd.read_excel('Supplementary Data 12.xlsx',sheet_name='function_annotation')
# print(df)
# 将Dataframe的每一行存为字典，ID列为key，average列为value
# df_dict = df.set_index('Geneid').to_dict()['Goid']
df_dict = df.set_index('#GeneID').to_dict()['Function']
# df_dict = df.set_index('#GeneID').to_dict()['GO_annotation']
# print(df_dict.keys())
# dd = pd.read_excel(r'C:\Users\dell\Desktop\Goldfish\data\flower_venn\显性共有.xlsx',sheet_name='only_gene')
dd = pd.read_excel('Supplementary Data 12.xlsx',sheet_name='GO_list')
# print(dd)

### 以下是KEGG步骤
# dff =pd.read_excel('Supplementary Data 12.xlsx',sheet_name='KEGG.Ko')
# dff =pd.read_excel('Supplementary Data 12.xlsx',sheet_name='KEGG pathway')
# # dff_dict = dff.set_index('Gene_id').to_dict()['KO']
# dff_dict = dff.set_index('KO').to_dict()['#pathway']
# print(dff_dict)
# ddd = pd.read_excel(r'gene_keggko.xlsx')

for i in dd['Geneid']:
    if i in df_dict.keys():
        print(i,df_dict[i])
    else:
        print(i)

# print('='*60)

# for j in dd['homo']:
#     if j in df_dict.keys():
#         print(j,df_dict[j])
#     else:
#         print(j)

# for i in ddd['KO']:
#     if i in dff_dict.keys():
#         print(i,dff_dict[i])
#     else:
#         print(i)