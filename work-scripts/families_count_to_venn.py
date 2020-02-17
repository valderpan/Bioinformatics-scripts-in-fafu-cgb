# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

# 把orthomcl生成的groups文件做venn图
import os
import re

os.chdir(r'D:\Result\师兄滴！\师兄做的orthomcl')
A='At'
C='Cs'
L='Longan'
V='Vv'
X='Xs'

# At = open('At.txt','w')
# Cs = open('Cs.txt','w')
# Longan = open('Longan.txt','w')
# Vv = open('Vv.txt','w')
# Xs = open('Xs.txt','w')

with open('groups.txt') as f :
    for line in f:
        line = line.strip()
        # if A in line:
        #     ID = re.findall('(OR_[0-9]+)',line)
        #     print('\n'.join(ID))

        # if C in line:
        #     ID = re.findall('(OR_[0-9]+)',line)
        #     print('\n'.join(ID))

        if L in line:
            ID = re.findall('(OR_[0-9]+)',line)
            print('\n'.join(ID))

        # if V in line:
        #     ID = re.findall('(OR_[0-9]+)',line)
        #     print('\n'.join(ID))

        # if X in line:
        #     ID = re.findall('(OR_[0-9]+)',line)
        #     print('\n'.join(ID))
