# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

'''
读取fasta文件并将其设置为不换行显示
计算碱基长度，计算'A''T''G''C'各自出现的次数！
'''
chr_gene = ''
with open(r'D:\Result\chr1.fasta','r') as f:
    for line in f:
        if line[0] != '>':
            line = line.strip()
            print(line)
        chr_gene=chr_gene +line
print(len(chr_gene))
A=chr_gene.count('A')
G=chr_gene.count('G')
C=chr_gene.count('C')
T=chr_gene.count('T')
print('the A counts is {}'.format(A))
print('the G counts is {}'.format(G))
print('the C counts is {}'.format(C))
print('the T counts is {}'.format(T))
