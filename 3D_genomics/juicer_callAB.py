#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/9/24


import sys
import argparse
sys.path.append('/share/home/stu_panhaoran/scripts')
import richlog as rl
import runshell as rs
import pandas as pd
from rich.traceback import install

__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20210924'
__version__ = 'v1.0'


def make_windows(genome,window):
    rs.run_command_noprint('bedtools makewindows -g {} -w {} > {}.windows'.format(genome,window,genome))


def read_ID(ChrID):
    ID = []
    with open(ChrID,'r') as f:
        lines = (line.strip() for line in f)
        for line in lines:
            ID.append(line)
    return ID


def run_jucier(ChrID,hicfile,window,juicer_path):
    for id in ChrID:
        rs.run_command('java -jar {} eigenvector KR {} {} BP {} {}.AB.txt'.format(juicer_path,hicfile,id,window,id))


def read_AB(ABresult):
    score = []
    with open(ABresult,'r') as f:
        lines = (line.strip() for line in f)
        for line in lines:
            score.append(line)
    return score


def read_windows(genome):
    df = pd.read_table('{}.windows'.format(genome),sep='\t',names=['seqid','start','end'])
    return df


def matchbinAB(ID,window_df):
    for id in ID:
        query_df = window_df[window_df['seqid'] == id]
        query_score = read_AB('{}.AB.txt'.format(id))
        query_df['score'] = query_score
        query_df.to_csv('{}.AB.result',sep='\t',header=False,index=False)


def main(args):
    install()
    genome = args.genome
    window = args.window
    juicer_path = args.juicer_path
    hic = args.hic
    chrid = args.chr
    rl.info_out('Make windows...,Window size:[{}]'.format(window))
    make_windows(genome,window)
    rl.info_out('Read ChrID...')
    ID = read_ID(chrid)
    rl.info_out('Run juicer call compartment...')
    run_jucier(ID,hic, window, juicer_path)
    window_df = read_windows(genome)
    rl.info_out('Match bin and ABscore, then output result...')
    matchbinAB(ID,window_df)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Use juicer to call compartment and output to tidy data format''',
        usage="python {} -g genome.fasta -w window_size -j juicer.jar -h .hic -c ChrID.txt ".format(sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,
                                                                            __mail__, __date__, __version__))
    parser.add_argument('-g', '--genome', required=True, help='Input Genome fasta file')
    parser.add_argument('-w', '--window', required=True, type=int,help='Input Window Size')
    parser.add_argument('-j', '--juicer_path', required = True, help='Input the path of juicer tools')
    parser.add_argument('-i', '--hic', required = True, help='Input the .hic file')
    parser.add_argument('-c', '--chr', required = True, help='Input the ChrID file')
    args = parser.parse_args()
    main(args)