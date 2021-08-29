#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/8/24


import sys
import argparse
import filter_fasta as ff
from rich.traceback import install

__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20210824'
__version__ = 'v1.0'


def read_seqID(file):
    ID = []
    with open(file,'r') as f:
        lines = (line.strip() for line in f)
        for line in lines:
            ID.append(line)
    return ID


def ExtractandOut_seq(fastaD,ID,output):
    with open(output,'w') as w:
        for key in fastaD.keys():
            if key in ID:
                w.write(">{}\n".format(key))
                w.write(str(fastaD[key])+'\n')


def main(args):
    install()
    fasta = args.fasta
    IDfile = args.queryid
    out = args.output
    fastaD =ff.store_fasta(fasta)
    ID = read_seqID(IDfile)
    ExtractandOut_seq(fastaD,ID,out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Extract query seq by header''',
        usage="python {} -f fasta -m mode -k keywords -o output.fasta ".format(sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,
                                                                            __mail__, __date__, __version__))
    parser.add_argument('-f', '--fasta', required=True, help='Input fasta')
    parser.add_argument('-q', '--queryid', required=True, help='Input query seqid file')
    parser.add_argument('-o', '--output', required=True, help='Output the extracted fasta')
    args = parser.parse_args()
    main(args)