# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-


##将Ptri_protein中的蛋白序列合成一行显示，并去掉每套序列末尾的'*'
import os
import re

os.chdir(r'D:\Result\speciestree')

print('>Ptrichocarpa')
with open('Ptri_protein.txt') as f:
    for line in f:
        line = line.replace('\n', '')

        line = re.findall(r'[^\*]',line)

        result =(''.join(line))

        print(result,end='')

        # result = result.replace('\n','')
        # print(result,end='')
