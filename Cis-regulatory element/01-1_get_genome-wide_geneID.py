# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

#Output the fa file as a line of header and line of seq
import os
import sys

fa_dict = {}
with open(sys.argv[1]) as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            fa_key = line[1:]
            fa_dict[fa_key] = ''
        else:
            line = line.upper()
            fa_dict[fa_key] += line.replace('\n','')
for i in fa_dict.keys():
    print(i)

#step2
#for i in fa_dict.keys():
#    print(len(fa_dict[i]))   
