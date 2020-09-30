# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

# 从BLAST结果中提取找到的所有的基因ID
# python gain_blast_result.py LA_blast_monoploidy_result.txt LA_monoploidy_result
import os
import sys

os.chdir(r'D:\Result\Transit\blast找甘蔗基因')
output_result = open(sys.argv[2],'w')
with open(sys.argv[1]) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('Query='):
            output_result.write('>>>>>'+line[0:18]+'\n')
        if line.startswith('Soffic.'):
            # print(line[])
            line = line.split('  ')
            output_result.write(line[0]+'\n')