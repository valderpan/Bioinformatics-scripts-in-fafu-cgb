# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

import os
import pandas as pd
import numpy as np
import re
import csv

tabfile_path = '/public1/home/stu_panhaoran/CRE_biotask/v0806/split_1M/plantCARE_output/tab'
output_path = '/public1/home/stu_panhaoran/CRE_biotask/v0806/split_1M/plantCARE_output/onlygroupby'
os.chdir(tabfile_path)
tab_file_list = os.listdir()
#print(main_file_list)

for i in tab_file_list:
    #sub_tab_file = os.listdir()
    #print(sub_tab_file)
    # print(i)
    sub_output_name = re.findall('\_([0-9]+)\.',i)[0]
    #print(sub_output_name)
    df = pd.read_table(i, sep='\t', header=None,quoting=csv.QUOTE_NONE)
    df.dropna(axis=0, how='any', inplace=True)
    df.rename(
        columns={0: 'GeneID', 1: 'Site Name', 2: 'Sequence', 3: 'Position', 4: 'Length', 5: 'Strand', 6: 'Organism',
                 7: 'Function'}, inplace=True)
    #df = df[~df['Function'].str.contains('?')]
    df = df[ ~ (df.iloc[:,7]=='?')]
    df = df.groupby(['GeneID', 'Site Name', 'Sequence', 'Strand', 'Length', 'Organism', 'Function'])['Position'].apply(
        list)
    df = df.reset_index()
    #df.to_excel(output_path+'/'+sub_output_name+'.xlsx',header=True,index=None)
    df.to_csv(output_path+'/'+sub_output_name+'.tab',sep='\t',header=True,index=False)
