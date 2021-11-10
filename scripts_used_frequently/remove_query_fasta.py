#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/11/10


import sys
import argparse
import filter_fasta as ff
from rich import print
from rich.traceback import install

__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20211110'
__version__ = 'v1.0'


def read_seqID(file):
    ID = []
    with open(file,'r') as f:
        lines = (line.strip() for line in f)
        for line in lines:
            ID.append(line)
    print('{} genes are obtained'.format(len(ID)))
    return ID



def remove_query_seq(fastaD,ID):
    removedD = {}
    match = 0
    for key in fastaD.keys():
        if key not in ID:
            removedD[key] = fastaD[key]
            match +=1
        else:
            continue
    return removedD,match


def Outseq(removedD,output):
    with open(output,'w') as w:
        for key in removedD.keys():
            w.write(">{}\n".format(key))
            w.write(str(removedD[key]) + '\n')


def main(args):
    install()
    fasta = args.fasta
    IDfile = args.queryid
    out = args.output

    fastaD = ff.store_fasta(fasta)
    ID = read_seqID(IDfile)
    removedD,match = remove_query_seq(fastaD,ID)
    print('{} seqs are are stored in the dictionary'.format(len(fastaD)))
    print('{} seqs are remain'.format(len(removedD)))

    if len(removedD) == len(fastaD) - len(ID):
        Outseq(removedD,out)
    else:
        print('The number of seqs before and after deletion is not equal!!!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Remove query seq by header''',
        usage="python {} -f fasta -q queryid.txt -o output.fasta ".format(sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,
                                                                            __mail__, __date__, __version__))
    parser.add_argument('-f', '--fasta', required=True, help='Input fasta')
    parser.add_argument('-q', '--queryid', required=True, help='Input need remove seqid file')
    parser.add_argument('-o', '--output', required=True, help='Output the extracted fasta')
    args = parser.parse_args()
    main(args)