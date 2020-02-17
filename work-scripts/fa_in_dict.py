# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# 将fa存入字典 并对每个header修改后输出原fa格式！
import sys
import os
fa_dict = {}
key_list = ''
os.chdir(r'./')
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('>'):
            line = line.rstrip()
            fa_key = line[0:20]
            key_list += fa_key

            fa_dict[fa_key] = ''
        else:
            fa_dict[fa_key] +=line
for i in fa_dict.keys():
    print(i.rstrip())
    print(fa_dict[i].rstrip())

