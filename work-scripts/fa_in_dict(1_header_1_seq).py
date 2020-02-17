# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

#将fasta文件存入字典，并将其显示为一行header一行seq的形式
import os
import sys
os.chdir(r'D:\Result')

fa_dict = {}
with open(sys.argv[1]) as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            fa_key = line[0:]
            fa_dict[fa_key] = ''
        else:
            line = line.upper()
            fa_dict[fa_key] += line.replace('\n','')
#将fasta文件显示为一行header一行seq！！！！！！！！！
for i in fa_dict.keys():
    print(i)
    print(fa_dict[i])