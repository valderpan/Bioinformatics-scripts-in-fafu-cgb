#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/8/17


import re
import sys
import argparse
import richlog as rl
from Bio import SeqIO
from rich.traceback import install


__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20210817'
__version__ = 'v1.0'


def store_fasta(fasta):
    D = {}
    for seq in SeqIO.parse(fasta,'fasta'):
        D[seq.id] = seq.seq
    return D


def filter_start(D,keywords,output):
    with open(output,'w') as w1:
        for key in D.keys():
            if key.startswith(keywords):
                w1.write('>'+key+'\n')
                w1.write(str(D[key])+'\n')


def filter_end(D,keywords,output):
    with open(output,'w') as w2:
        for key in D.keys():
            if key.endswith(keywords):
                w2.write('>'+key+'\n')
                w2.write(str(D[key])+'\n')


def filter_re(D,keywords,output):
    with open(output,'w') as w3:
        for key in D.keys():
            kw = re.findall(keywords, key)[0]
            if kw in key:
                w3.write('>'+key+'\n')
                w3.write(str(D[key])+'\n')


def main(args):
    install()
    fasta = args.fasta
    kw = args.keywords
    mode = args.mode
    output = args.output
    rl.info_out('Read fasta and save to dictionary...')
    D = store_fasta(fasta)
    rl.info_out('Filtering fasta...')
    if mode == 'start':
        filter_start(D,kw,output)
    elif mode == 'end':
        filter_end(D,kw,output)
    elif mode == 're':
        filter_re(D,kw,output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Filter fasta by header''',
        usage="python {} -f fasta -m mode -k keywords -o output.fasta ".format(sys.argv[0]),
        epilog= 'author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,
															__mail__,__date__,__version__))
    parser.add_argument('-f','--fasta',required=True,help='Input fasta')
    parser.add_argument('-m','--mode',required=True,choices=['start','end','re'],help='Choose filter mode')
    parser.add_argument('-k','--keywords',required=True,type=str,help='Input the filtered keywords or Regular expressions')
    parser.add_argument('-o','--output',required=True,help='Output the filtered fasta')
    args = parser.parse_args()
    main(args)

