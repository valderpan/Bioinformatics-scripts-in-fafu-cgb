# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

# 将cds文件中名字里的‘ ’变为‘|’用于下一步操作！
import os

os.chdir(r'D:\databases\yellowhorn')
with open('Vvinifera_145_cds.fa') as f :
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if len(line) < 29:
                line =line
                # print(line)
            else:
                line = line[0:13]+'|'+line[14:]
        if not line.startswith('>'):
            line =line.upper()
        print(line)