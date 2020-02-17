# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

# 将Trinity组装之后的转录本翻译为蛋白质后，修改pep文件的header
# 使用：linux系统内 python * A.pep.fasta A.fasta
import os
import sys

os.chdir(r'D:\Result\Transit\转录组组装')
fileList = os.listdir('./')
# output = open(sys.argv[2],'w')
fa_dict = {}
fasta_number = 0

with open('A.pep.fasta') as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            fasta_number += 1
            fa_key = line[0:30] #若取值很短会造成键重复，从而造成丢失数据
            fa_dict[fa_key] = ''
        else:
            fa_dict[fa_key] += line
# print(fasta_number)

number_list = []
for i in range(1, fasta_number+1):
    i = '>' + 'A' + '%05d' % i
    number_list.append(i)
# print(number_list)

fa_dict_values = fa_dict.values()
# print(len(fa_dict_values))
new_fa_dict = dict(zip(number_list,fa_dict_values))

for j in new_fa_dict.keys():
    print(j)
    print(new_fa_dict[j])

# output.close()