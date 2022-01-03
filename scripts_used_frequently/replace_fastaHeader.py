#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2022/1/3


'''
%prog <oldID.fasta> <IDlink.txt> > <newID.txt>

Replace the header of fasta
Note: IDlink.txt must be a file split by "\\t"

>>> python %prog .pep.fasta gene_rename.txt(split by "\\t") > .pep.newID.fasta

__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20220103'
__version__ = 'v1.0'

'''


import sys
import fa2dict_Bio as fdB


def read_replacefile(file):
    old2new = {}
    with open(file) as f:
        lines = (line.strip() for line in f)
        for line in lines:
            oldID,newID = line.split('\t')
            old2new[oldID] = newID
    return old2new


def Replace_fastaHeader(old2new,fastaD):
    ReplaceD = {}
    for key in fastaD.keys():
        ReplaceD[old2new[key]] = fastaD[key]
    return ReplaceD


def Sort_output_Fasta(fastaD):
    for key in sorted(fastaD.items(),key=lambda x:x[0]):
        print(">" + key[0])
        print(str(fastaD[key[0]]))


if __name__ == '__main__':
    from optparse import OptionParser
    from rich.traceback import install
    install()
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args) == 2:
        fasta_file = args[0]
        IDlink = args[1]
        fastaD = fdB.Fa2dict_Bio(fasta_file)
        old2new = read_replacefile(IDlink)
        Replace_D = Replace_fastaHeader(old2new,fastaD)
        Sort_output_Fasta(Replace_D)
    else:
        sys.exit(p.print_help())
