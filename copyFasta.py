#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/09/23

'''
%prog <original fasta file> <ctgID file> <output file>

Copy the corresponding contig fasta infomation
After collapse_rescue, a new contig with the _d suffix is generated, 
and a copy of the new contig is made according to the previous contig fasta information

>>> python %prog original.fasta ctgID.txt output

__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20220924'
__version__ = 'v1.0'

'''

import sys
from Bio import SeqIO
from rich.traceback import install

def Fa2dict(fasta_file):
    fastaD = {}
    for seq in SeqIO.parse(fasta_file,'fasta'):
        fastaD[seq.id] = seq.seq
    return fastaD

def copyfasta(fastaD,ctgID,output):
    ctgL = []
    with open(ctgID) as f:
        lines = (line.strip() for line in f)
        for line in lines:
            ctgL.append(line)

    with open(output,'w') as w1:
        for i in ctgL:
            if i in fastaD.keys():
                w1.write('>'+i+'\n')
                w1.write(str(fastaD[i])+'\n')
            else:
                newi = i.split('_')[0]
                # print(newi)
                w1.write('>'+i+'\n')
                w1.write(str(fastaD[newi])+'\n')


if __name__ == "__main__":
    from optparse import OptionParser
    from rich.traceback import install
    install()
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args)==3:
        fasta_file = args[0]
        ctg_file = args[1]
        output = args[2]
        fastaD= Fa2dict(fasta_file)
        copyfasta(fastaD,ctg_file,output)
    else:
        sys.exit(p.print_help())