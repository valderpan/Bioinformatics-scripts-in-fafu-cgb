# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

import os
os.chdir(r'D:\Result\test')
codon_table = {
    'GCU':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'CGU':'R', 'CGC':'R',
    'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R', 'UCU':'S', 'UCC':'S',
    'UCA':'S', 'UCG':'S', 'AGU':'S', 'AGC':'S', 'AUU':'I', 'AUC':'I',
    'AUA':'I', 'UUA':'L', 'CUU':'L', 'CUC':'L', 'CUA':'L',
    'CUG':'L', 'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'GUU':'V',
    'GUC':'V', 'GUA':'V', 'ACU':'T', 'ACC':'T', 'ACA':'T',
    'ACG':'T', 'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'AAU':'N',
    'AAC':'N', 'GAU':'D', 'GAC':'D', 'UGU':'C', 'UGC':'C', 'CAA':'Q',
    'CAG':'Q', 'GAA':'E', 'GAG':'E', 'CAU':'H', 'CAC':'H', 'AAA':'K',
    'AAG':'K', 'UUU':'F', 'UUC':'F', 'UAU':'Y', 'UAC':'Y',
    'UGG':'W',
    'UAG':'STOP', 'UGA':'STOP', 'UAA':'STOP','AUG':'START','GUG':'START','UUG':'START'
    }

rna = ''
with open('rna.txt') as f :
    for line in f :
        if not line.startswith('>'):
            rna = rna + line.strip()
    print(len(rna))
    for frame in range(3):
        prot = ''
        print('reading frame:'+str(frame+1))
        print(frame)
    for i in range(0,len(rna),3): # 这里若将0换为frame，则是从2-186，每隔3个输出一次
        conda=rna[i:i+3]
        print(conda)
        if conda in codon_table:
            if codon_table[conda] == 'STOP':
                prot = prot + '*'
            else:
                prot = prot + codon_table[conda]
        else:
            prot = prot+'-'
    print(prot)

        # ii=0
        # while ii < 48 :
        #     print(prot[i:i+48])
        #     i = i +48
