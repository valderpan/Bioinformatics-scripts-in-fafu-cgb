# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

# 从文件中找到特定位置的基因用作blast文件进行比对
import os

fa_D = {}
os.chdir(r'D:\Result\hg19')
with open('Sspon.HiC_chr_asm.fasta') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            chr_name = line[1:]
            fa_D[chr_name] = ''
        else:
            fa_D[chr_name] = line


print(fa_D['Chr7D'][12144100:12146398])
print(fa_D['Chr7A'][13631480:13633765])
print(fa_D['Chr7A'][13646544:13648835])
print(fa_D['Chr7C'][9013327:9015613])
print(fa_D['Chr7C'][8940410:8941030])
print(fa_D['Chr7C'][12549487:12551184])
print(fa_D['Chr7B'][12559924:12561111])
print(fa_D['Chr7B'][12551425:12552011])



