#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan

from Bio import SeqIO

def Fa2dict_Bio(fasta_file):
    fasta_dict = {}
    for seq in SeqIO.parse(fasta_file,'fasta'):
        fasta_dict[seq.id] = seq.seq
    return fasta_dict
