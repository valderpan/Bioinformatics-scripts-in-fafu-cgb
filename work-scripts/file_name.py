# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# 将一个目录下的文件名提取出来(.list前面的)
import os
import re
os.chdir(r'D:\Result\test\test')
file = os.listdir()
for i in file:
    if  i[-4:] == 'list':
        name = re.findall('([a-zA-Z0-9]+)\.',i)
        for file_name in name:
            print(file_name)