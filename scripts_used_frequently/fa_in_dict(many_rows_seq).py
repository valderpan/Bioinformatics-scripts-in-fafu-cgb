# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# 将一个fasta文件保存为字典，且保留其原格式！
# 当一个fasta文件中的seq有多行时使用这个脚本！
import os
import sys
from collections import OrderedDict
Dict = OrderedDict()
# os.chdir(r'D:\Result')
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('>'):
            line = line.rstrip()
            Dict_key = line[1:]
            Dict[Dict_key] = ''
        else:
            Dict[Dict_key] += line
    print(Dict['S_chr14_97321636_17M824N81M309N32M/1'])