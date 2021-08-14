#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/8/12

'''
统计NpX与AP85-441之间共线区域的AB swith、保守情况
'''

import re
import sys
import gff2bed
import argparse
import pandas as pd
from path import Path
from rich.console import Console
from rich.traceback import install
from intervaltree import Interval,IntervalTree


__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20210813'
__version__ = 'v1.0'

def concatABfile(pcafile_path):
    files = [i for i in Path(pcafile_path).files() if i.endswith('.compartment.csv')]
    L = []
    for file in files:
        df = pd.read_csv(file)
        L.append(df)
    c_df = pd.concat(L)
    return c_df


def storeChrinInter(c_df):
    tree_Chr = {}
    for index,row in c_df.iterrows():
        if row['score'] < 0:
            AB = 'B'
        elif row['score'] > 0:
            AB = 'A'
        else:
            AB = 'NaN'

        if row['seqnames'] not in tree_Chr:
            tree_Chr[row['seqnames']] = IntervalTree()
            tree_Chr[row['seqnames']].addi(row['start'],row['end'],AB)
        else:
            tree_Chr[row['seqnames']].addi(row['start'], row['end'], AB)
    return tree_Chr


def DictchrInter(c_df):
    Chr_D = {}
    for seq in c_df['seqnames'].unique().tolist():
        seq_D = {}
        query_df = c_df[c_df['seqnames']==seq]
        for index,row in query_df.iterrows():
            if row['seqnames'] not in Chr_D:
                Chr_D[row['seqnames']] = {}
                seq_D[row['end']]= [row['start'],row['end']]
                Chr_D[row['seqnames']].update(seq_D)
            else:
                seq_D[row['end']] = [row['start'], row['end']]
                Chr_D[row['seqnames']].update(seq_D)
    return Chr_D

def store_gene(g_bed):
    tree_gene = {}
    for index,row in g_bed.iterrows():
        chr,start,end,ID = row
        if chr not in tree_gene:
            tree_gene[chr] = IntervalTree()
            tree_gene[chr].addi(start,end,ID)
        else:
            tree_gene[chr].addi(start,end,ID)
    return tree_gene


def classifyGene2AB(tree_Chr,tree_gene):
    for key in tree_gene.keys():
        for inter in sorted(tree_gene[key]):
            b_inter = tree_Chr[key].overlap(inter)
            for b in b_inter:
                if inter.begin >= b.begin and inter.end < b.end: #完全包含在内
                    n_inter = Interval(inter.begin, inter.end,'@'.join([inter.data,b.data]))
                    tree_gene[key].remove(inter)
                    tree_gene[key].add(n_inter)
                elif inter.end > b.begin and inter.begin < b.begin: #卡在左边的
                    if inter.end - b.begin >= (inter.end-inter.begin)/2: #卡在左边但是整体是位于本体内的
                        n_inter = Interval(inter.begin, inter.end, '@'.join([inter.data, b.data]))
                        tree_gene[key].remove(inter)
                        tree_gene[key].add(n_inter)
                elif inter.begin < b.end and inter.end > b.end:
                    if b.end -inter.begin > (inter.end - inter.begin)/2: #卡在右边但是整体是位于本体内的
                        n_inter = Interval(inter.begin, inter.end, '@'.join([inter.data, b.data]))
                        tree_gene[key].remove(inter)
                        tree_gene[key].add(n_inter)
    ced_t_gene = tree_gene
    return ced_t_gene


def Output(ced_t_gene,output):
    with open(output,'w') as w:
        for key in ced_t_gene.keys():
            for inter in sorted(ced_t_gene[key]):
                ID,AB = inter.data.split('@')
                w.write('{}\t{}\t{}\t{}\t{}\n'.format(key,inter.begin,inter.end,ID,AB))


def main(args):
    console = Console()
    gff = args.gff3
    Type = args.type
    ABpath = args.path
    output = args.output
    install()
    console.log('Read the gff3 file...')
    df, pre = gff2bed.read_gff(gff)
    console.log('Convert gff to bed and output result...')
    g_bed = gff2bed.Convertgff2bed(df, Type)
    gff2bed.output(g_bed,pre,Type)
    console.log('Remove contig information...')
    cl_bed = gff2bed.removectg(g_bed)
    console.log('Read ABcompartment files...')
    c_df = concatABfile(ABpath)

    Chr_D = DictchrInter(c_df)
    console.log('Constructing an intervaltree at the chromosome bin resolution level...')
    tree_Chr = storeChrinInter(c_df)
    console.log('Constructing an intervaltree of gene start and stop positions...')
    tree_gene = store_gene(cl_bed)
    console.log('Classify the target genes into ABcompartment...')
    ced_t_gene = classifyGene2AB(tree_Chr, tree_gene)
    Output(ced_t_gene,output)
    console.log('Done!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Classify each gene into AB compartment''',
        usage="python {} -g gff3 -t gffe_type -p path -o output.bed ".format(sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,
                                                                            __mail__, __date__, __version__))
    parser.add_argument('-g', '--gff3', required=True, help='Input whole genome gff3 file')
    parser.add_argument('-t', '--type', required=True,
                        help='Select the type of sequence to be extracted from the gff3 file',
                        choices=['gene','CDS','mRNA'])
    parser.add_argument('-p', '--path', required=True,
                        help='Input the path where the result file of the call compartment is located')
    parser.add_argument('-o', '--output', required=True,
                        help='Output the final result (bed format)')
    args = parser.parse_args()
    main(args)