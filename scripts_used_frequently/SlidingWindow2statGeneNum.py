#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan

import sys
import math
import argparse
from Bio import SeqIO
from intervaltree import IntervalTree, Interval

__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20210729'
__version__ = 'v1.0'


def read_genome(genome):
    geno_D = {}
    for seq in SeqIO.parse(genome, 'fasta'):
        geno_D[seq.id] = seq.seq
    return geno_D


def read_gff_intree(gff3):
    '''
    将每个gene的位置存入字典树
    :param gff3:
    :return:
    '''
    chr_tree = {}
    with open(gff3,'r') as f:
        lines = (line.strip() for line in f)
        for line in lines:
            line_list = line.split('\t')
            if line_list[2] == 'gene':
                geneID = line_list[8].split(';')[0].split('=')[1]
                if line_list[0] not in chr_tree:
                    chr_tree[line_list[0]] = IntervalTree()
                    chr_tree[line_list[0]].addi(int(line_list[3]),int(line_list[4]),geneID)
                else:
                    chr_tree[line_list[0]].addi(int(line_list[3]), int(line_list[4]), geneID)
    return chr_tree


def slid_window(geno_D,window,step):
    '''
    将每条序列按照window和step的大小进行划分和移动，并存为区间树
    :param geno_D:
    :param window:
    :param step:
    :return:
    '''
    tree_D = {}
    for key in geno_D.keys():
        windows = math.ceil(len(geno_D[key])/step)
        win_s = 0
        for i in range(windows):
            if i == 0 :
                win_e = win_s + window
                if key not in tree_D:
                    tree_D[key] = IntervalTree()
                    tree_D[key].addi(win_s,win_e,0)
                else:
                    tree_D[key].addi(win_s,win_e,0)
            else:
                if i < windows-1:
                    win_s += step
                    win_e = win_s + window

                    tree_D[key].addi(win_s,win_e,0)
                else:
                    win_s += step
                    win_e = len(geno_D[key])
                    tree_D[key].addi(win_s, win_e,0)

    return tree_D

def parse_overlap(chr_tree,tree_D):
    '''
    将序列区间树与gff位置树进行匹配
    :param chr_tree:
    :param tree_D:
    :return:
    '''
    for key in tree_D.keys():

        for inter in sorted(tree_D[key]):
            all_cont = 0
            if chr_tree[key].overlap(inter):

                for i in chr_tree[key].overlap(inter):
                    if i.begin >= inter.begin and i.end < inter.end:
                        all_cont +=1
                    elif i.end >= inter.begin and i.begin < inter.begin:
                        if i.end-inter.begin >= (i.end-i.begin)/2:
                            all_cont += 1

                        else:
                            continue
                    elif i.end > inter.end and i.begin <= inter.end:
                        if i.end-inter.end >= (i.end-i.begin)/2:
                            continue
                        else:
                            all_cont += 1
            new_inter = Interval(inter.begin,inter.end,inter.data+all_cont)
            tree_D[key].remove(inter)
            tree_D[key].add(new_inter)

    return tree_D


def output_result(newtree_D):
    for key in sorted(newtree_D.keys()):
        for inter in sorted(newtree_D[key]):
            print(key,inter.begin,inter.end,inter.data,sep='\t')


def main(args):
    fasta_file = args.fasta
    gff_file = args.gff
    window = args.window
    step = args.step

    geno_D = read_genome(fasta_file)
    chr_tree = read_gff_intree(gff_file)
    tree_D = slid_window(geno_D, window, step)
    newtree_D = parse_overlap(chr_tree, tree_D)
    output_result(newtree_D)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Sliding window to stat gene num on genome-wide fasta''',
        usage="python {} -f fasta -g gff -w window -s step > output.bed".format(sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__, __mail__, __date__,
                                                                            __version__))
    parser.add_argument('-f', '--fasta', required=True, help='Input fasta file')
    parser.add_argument('-g', '--gff', required=True, help='Input gff3 file')
    parser.add_argument('-w', '--window', required=True,type=int, help='Specify the window size')
    parser.add_argument('-s', '--step', required=True,type=int, help='Specify the step size')
    args = parser.parse_args()
    main(args)
