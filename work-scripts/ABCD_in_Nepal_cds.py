# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

# 将尼珀尔全cds序列中的ABCD等位基因的cds提取出来
import os
import sys
os.chdir(r'D:\Result\Transit\MCscan\Nipor')

fa_dict = {}
A_result = open('Nipor_Ss.v20191222.A_cds.fasta','w')
B_result = open('Nipor_Ss.v20191222.B_cds.fasta','w')
C_result = open('Nipor_Ss.v20191222.C_cds.fasta','w')
D_result = open('Nipor_Ss.v20191222.D_cds.fasta','w')

with open('Nipor_Ss.v20191222.cds.fasta') as f:
    for line in f:
        line =line.rstrip()
        if line.startswith('>'):
            fa_key = line[0:]
            fa_dict[fa_key] = ''
        else:
            fa_dict[fa_key] += line
#
for i in fa_dict.keys():
    if i[8] == 'A':
        A_result.write(i+'\n')
        A_result.write(fa_dict[i]+'\n')
    if i[8] == 'B':
        B_result.write(i+'\n')
        B_result.write(fa_dict[i]+'\n')
    if i[8] == 'C':
        C_result.write(i+'\n')
        C_result.write(fa_dict[i]+'\n')
    if i[8] == 'D':
        D_result.write(i+'\n')
        D_result.write(fa_dict[i]+'\n')
