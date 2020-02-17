# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

# 根据筛选出来的contigID对CDS进行过滤，过滤掉contig序列！
import os

Xs_dict = {}

os.chdir(r'D:\databases\yellowhorn')
#将Xs.cds中的序列存入字典！
with open('yellowhorn.cds') as X:
    for line in X:
        line =line.strip()
        if line.startswith('>'):
            Xs_key = line[1:]
            Xs_dict[Xs_key] = ''
        else:
            Xs_dict[Xs_key] += line.replace('\n','')
#进行筛选！
contig_id = []
with open('yellowhorn_contig_list.txt') as f :
    for i in f :
        i = i.strip()
        contig_id.append(i)
    for key in Xs_dict.keys():
        if not key in contig_id:
            print('>'+key)
            print(Xs_dict[key])
