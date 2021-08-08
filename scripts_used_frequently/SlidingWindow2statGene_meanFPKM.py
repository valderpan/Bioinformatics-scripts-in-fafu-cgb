#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan

import sys
import SlidingWindow2statGeneNum as SW
import SlidingWindow2statGeneID as SWI
import argparse
import pandas as pd


__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20210806'
__version__ = 'v1.0'


def store_gene2fpkm(fpkm_file):
    df = pd.read_excel(fpkm_file,sep='\t')
    c = [column for column in df]
    t_D = {}
    for i in c[1:]:
        tmp_D = df.set_index(c[0]).to_dict()[i]
        t_D[i] = tmp_D
    return t_D,df


def match_bed2fpkm(newtree_D,t_D,fpkm_df):
    for key in sorted(newtree_D.keys()):
        for inter in sorted(newtree_D[key]):
            fpkm_L = []
            for i in t_D.keys():
                q_df = fpkm_df[fpkm_df['new'].isin(inter.data)]
                mean_fpkm = round(q_df.loc[:,i].sum()/q_df.shape[0],3)
                fpkm_L.append(str(mean_fpkm))

            print(key, inter.begin, inter.end, '\t'.join(fpkm_L), sep='\t')



def main(args):
    fasta_file = args.fasta
    gff_file = args.gff
    window = args.window
    step = args.step
    fpkm_file = args.fpkm

    geno_D = SW.read_genome(fasta_file)
    chr_tree = SWI.read_gff_intree(gff_file)
    tree_D = SW.slid_window(geno_D, window, step)
    newtree_D = SWI.parse_overlap(chr_tree, tree_D)
    t_D,fpkm_df=store_gene2fpkm(fpkm_file)
    match_bed2fpkm(newtree_D,t_D,fpkm_df)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Sliding window to stat gene mean fpkm''',
        usage="python {} -f fasta -g gff -w window -s step > output.bed".format(sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__, __mail__, __date__,
                                                                            __version__))
    parser.add_argument('-f', '--fasta', required=True, help='Input fasta file')
    parser.add_argument('-g', '--gff', required=True, help='Input gff3 file')
    parser.add_argument('-p', '--fpkm', required=True, help='Input fpkm file')
    parser.add_argument('-w', '--window', required=True,type=int, help='Specify the window size')
    parser.add_argument('-s', '--step', required=True,type=int, help='Specify the step size')

    args = parser.parse_args()
    main(args)
