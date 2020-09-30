# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

'''
需求：查看每个bin里面有多少reads mapping到compartment中的
使用：python RnaSeq_bam_bins_readsmapping.py xxx.mapping.position1.txt xxx/xxx/AB_region/(绝对路径) > xxx.RnaSeq_bins_readsmappingstate.txt
'''

import pandas as pd
import os
import re
from interval import Interval
import math


def get_bam_info(bam_file,AB_region_file_path):
    df = pd.read_table(bam_file, sep='\t', header=None)
    df = df[df.iloc[:, 0].str.contains('Chr')]
    print(df)
    os.chdir(AB_region_file_path)
    file_list = os.listdir(AB_region_file_path)
    for file in file_list:
        file_chr_name = re.findall('([0-9a-zA-Z]+)\.',file)[0][:-1]    #TODO 注意这里
        # print(file_chr_name)
        df1 = df[df.iloc[:,0].str.contains(file_chr_name)]
        # print(df1)
        value_list = df1.iloc[:, 1].tolist()
        # print(value_list)
        AB_region_inter_list = []
        with open(file) as f:
            # print(file)
            for line in f:
                line = line.strip()
                if line.startswith('start'):
                    continue
                inter_start,inter_end = line.split()
                inter = Interval(int(inter_start), int(inter_end))
                AB_region_inter_list.append(inter)
        # print(AB_region_inter_list)
                num=0
                for i in value_list:
                    if i in inter:
                        num +=1
                if num != 0:
                    print(file_chr_name,num,sep='\t')

get_bam_info(r'C:\Users\dell\Desktop\脚本测试\20200819 RNA-Seq_bam_bins_readsa_mapping\test.mapping.position.txt',
             r'C:\Users\dell\Desktop\脚本测试\20200819 RNA-Seq_bam_bins_readsa_mapping\AB_region')