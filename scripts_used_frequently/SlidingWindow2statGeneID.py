#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan

import sys
import SlidingWindow2statGeneNum as SW
import argparse
from intervaltree import IntervalTree, Interval

__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20210806'
__version__ = 'v1.0'


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
                if line_list[0] not in chr_tree:
                    chr_tree[line_list[0]] = IntervalTree()
                    chr_tree[line_list[0]].addi(int(line_list[3]),int(line_list[4]),[line_list[8].split(';')[0].split('=')[1]])
                else:
                    chr_tree[line_list[0]].addi(int(line_list[3]), int(line_list[4]), [line_list[8].split(';')[0].split('=')[1]])
    return chr_tree


def parse_overlap(chr_tree,tree_D):
    '''
    将序列区间树与gff位置树进行匹配
    :param chr_tree:
    :param tree_D:
    :return:
    '''
    for key in tree_D.keys():

        for inter in sorted(tree_D[key]):
            all_cont = []
            if chr_tree[key].overlap(inter):

                for i in chr_tree[key].overlap(inter):
                    if i.begin >= inter.begin and i.end < inter.end:
                        all_cont.extend(i.data)
                    elif i.end >= inter.begin and i.begin < inter.begin:
                        if i.end-inter.begin >= (i.end-i.begin)/2:
                            all_cont.extend(i.data)
                        else:
                            continue
                    elif i.end > inter.end and i.begin <= inter.end:
                        if i.end-inter.end >= (i.end-i.begin)/2:
                            continue
                        else:
                            all_cont.extend(i.data)
            new_inter = Interval(inter.begin,inter.end,all_cont)
            tree_D[key].remove(inter)
            tree_D[key].add(new_inter)

    return tree_D

def output_result(newtree_D):
    for key in sorted(newtree_D.keys()):
        for inter in sorted(newtree_D[key]):
            print(key,inter.begin,inter.end,','.join(inter.data),sep='\t')


def main(args):
    fasta_file = args.fasta
    gff_file = args.gff
    window = args.window
    step = args.step

    geno_D = SW.read_genome(fasta_file)
    chr_tree = read_gff_intree(gff_file)
    tree_D = SW.slid_window(geno_D, window, step)
    newtree_D = parse_overlap(chr_tree, tree_D)
    output_result(newtree_D)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Sliding window to stat gene ID''',
        usage="python {} -f fasta -g gff -w window -s step > output.bed".format(sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__, __mail__, __date__,
                                                                            __version__))
    parser.add_argument('-f', '--fasta', required=True, help='Input fasta file')
    parser.add_argument('-g', '--gff', required=True, help='Input gff3 file')
    parser.add_argument('-w', '--window', required=True,type=int, help='Specify the window size')
    parser.add_argument('-s', '--step', required=True,type=int, help='Specify the step size')
    args = parser.parse_args()
    main(args)