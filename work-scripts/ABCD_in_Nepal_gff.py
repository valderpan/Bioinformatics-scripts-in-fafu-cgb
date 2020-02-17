# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

# 将尼珀尔全GFF序列中的ABCD等位基因的gff提取出来
import os
import pandas as pd

os.chdir(r'D:\Result\Transit\MCscan\Nipor')

A_result = open('Nipor_Ss.v20191222.A_Chr.gff3','w')
B_result = open('Nipor_Ss.v20191222.B_Chr.gff3','w')
C_result = open('Nipor_Ss.v20191222.C_Chr.gff3','w')
D_result = open('Nipor_Ss.v20191222.D_Chr.gff3','w')

with open('Nipor_Ss.v20191222.Chr.gff3') as f:
    for line in f:
        line =line.rstrip()
        if line[5] == 'A' or line[4] == 'A':
            A_result.write(line+'\n')
        if line[5] == 'B' or line[4] == 'B':
            B_result.write(line+'\n')
        if line[5] == 'C' or line[4] == 'C':
            C_result.write(line+'\n')
        if line[5] == 'D' or line[4] == 'D':
            D_result.write(line+'\n')