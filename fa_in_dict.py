# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

## 将fa文件保存为dict，其中‘>’开头的为key，seq为value
import os
os.chdir(r'D:\Result')
genome = {}
read_num=[]
read_seq=[]
with open('reads_1.fa') as input_fa:
    for line in input_fa:
        if line[0] == '>':
            line =line.strip()
            read_num.append(line[1:])
        else:
            line =line.strip()
            read_seq.append(line)

genome=dict(zip(read_num,read_seq))

print(genome['2'][0:10])
print(genome['777'][0:10])
