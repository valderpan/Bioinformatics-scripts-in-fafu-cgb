#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/03/06


import sys
import fa2dict_Bio as fdB
import richlog as rl
from rich import print

__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20220306'
__version__ = 'v1.0'


def stat_seqLen(fastaD):
    fasta2Len = {}
    for key in fastaD.keys():
        fasta2Len[key] = len(fastaD[key])
        print('{}\t{}'.format(key,len(fastaD[key])))
    return fasta2Len


def sliding_windows(fastaD,window,step,output):
    with open(output,'w') as w:
        for key in fastaD.keys():
            start = 1
            while start < fastaD[key]:
                if start+window < fastaD[key]:
                    w.write('{}\t{}\t{}\n'.format(key,start,start+window-1))
                    start += step
                else:
                    w.write('{}\t{}\t{}\n'.format(key,start,fastaD[key]))
                    break


def main(args):
    fasta = args.fasta
    window = args.window
    step = args.step
    output = args.output

    rl.info_out('Step1: Count the length of each sequence in the fasta file...')
    fasta_Dict = fdB.Fa2dict_Bio(fasta)
    fasta2Len  = stat_seqLen(fasta_Dict)
    rl.info_out('Step2: Create window files...')
    sliding_windows(fasta2Len,window,step,output)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''This script is used to create the window file''',
        usage="python {} -f genome.fasta -w window_size -s step_size -o output ".format(sys.argv[0]),
        epilog= 'author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
    parser.add_argument('-f','--fasta',required=True,help='Input fasta file')
    parser.add_argument('-w','--window',required=True,type=int,help='Specify window size')
    parser.add_argument('-s','--step',required=True,type=int,help='Specify step length')
    parser.add_argument('-o','--output',required=True,help='Specify the output file name')
    args = parser.parse_args()
    main(args) 