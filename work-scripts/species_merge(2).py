# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

import os
#import re
os.chdir(r'D:\Result\todomafft\Athaliana')

# fw = open('out.fasta','w')
Ptri_genome = {}

##读取Ptri的全部的蛋白序列，并将其保存为字典
with open('Athaliana_447_Araport11.protein_primaryTranscriptOnly.fa') as Ptri_fa:

    for line in Ptri_fa:
         if line[0] == '>':
             name=line.split()[0][1:]
             Ptri_genome[name]=''
         else:
            Ptri_genome[name] += line.replace('\n','')
# print(Ptri_genome)
# print(Ptri_genome.keys())
# print(Ptri_genome.values())

with open('AthaID.txt') as f:
    for i in f:
        i = i.rstrip()
        print(i)
        for key in Ptri_genome.keys():
            if i == key:
                print(Ptri_genome[i])
