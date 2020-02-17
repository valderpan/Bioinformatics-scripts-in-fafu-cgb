# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-


import os
#import re
os.chdir(r'D:\Result\speciestree')

# fw = open('out.fasta','w')
# Ptri_genome = {}
#
# ##读取Ptri的蛋白序列，并将其保存为字典
# with open('Ptrichocarpa_444_v3.1.protein_primaryTranscriptOnly.fa') as Ptri_fa:
#
#     for line in Ptri_fa:
#          if line[0] == '>':
#              name=line.split()[0][1:18]
#              Ptri_genome[name]=''
#          else:
#             Ptri_genome[name] += line.replace('\n','')

    # print(Ptri_genome)
    # print(Ptri_genome.keys())
    # for i in Ptri_genome.keys():
    #     print(i)
    # print('-'*50)

with open('species_merge_Ptri.txt') as f:

    for line in f:
        line = line.strip()

        # print(line)
        lines = line.split(',')
        # print(lines)
        for num in lines :
            # print(num)
            if num.startswith(' '):
                num = num[1:]
            if num.startswith('Ptri'):
                num = num[13:]
            print(num)

## 得到单个物种的单拷贝基因ID
## 下一步根据单拷贝基因ID在蛋白文件中将对应的蛋白序列提出来！


        # print(i)
        # for i in Ptri_genome.keys():
        #     print(i)
        #     if i == num:
        #         print(i)
        # print('-'*50)
# for i in Ptri_genome.keys():
#     # if i.startswith('Potri.T155200.1.p'):
#     # print(num)
#     print(i)
#     print(Ptri_genome[i])
#     if i == num:
#         print(i)
#             fw.write(i)
#             fw.write('\n')
#             fw.write(Ptri_genome[i])
#             fw.write('\n')
# fw.close()


# for i in Ptri_genome.keys():
#     print(i)