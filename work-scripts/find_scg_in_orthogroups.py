# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

import os

gene_dict ={}
os.chdir(r'D:\databases\test\Orthogroups')
with open('Orthogroups.txt') as  f:
    for line in f:
        line = line.strip()
        # gene_dict[name] =
        if not line.startswith('Orthogoup'):
            id = line[0:9]
            # print(line[11:])
            gene_dict[id] = line[11:]
# print(gene_dict)

with open('Orthogroups_SingleCopyOrthologues.txt') as l:
    for i in l:
        i = i.strip()
        print(i,gene_dict[i])