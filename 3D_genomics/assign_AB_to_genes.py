#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/12/2


__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20211203'
__version__ = 'v1.0'


import sys
import pandas as pd
from intervaltree import Interval,IntervalTree


def read_compartmentfile(file):
    comdf = pd.read_csv(file)
    return comdf


def read_anchors(file):
    df = pd.read_table(file,sep='\t',comment='#',names=['s1id','s2id','score'])
    s1id = df['s1id'].unique().tolist()
    s2id = df['s2id'].unique().tolist()
    return s1id,s2id


def read_allgenebed(file):
    df = pd.read_table(file,sep='\t',names=['seqid','start','end','ID','uk1','uk2'])
    return df


def Find_query_genebed(s1id,s2id,s1genebed,s2genebed):
    # print(s1id[1])
    # print(s2id[1])
    # print(s1genebed)
    # print(s2genebed)
    # print('='*50)
    s1qdf = s1genebed[s1genebed['ID'].isin(s1id)]
    s2qdf = s2genebed[s2genebed['ID'].isin(s2id)]
    # print(s1qdf)
    # print(s2qdf)
    if len(s1id)== s1qdf.shape[0] and len(s2id) == s2qdf.shape[0]:
        print('All corresponding genes have been extracted ~')
    else:
        print('s1id len : {}'.format(len(s1id)))
        print('s2id len : {}'.format(len(s2id)))
        print('s1qdf len : {}'.format(s1qdf.shape[0]))
        print('s2qdf len : {}'.format(s2qdf.shape[0]))
        print('Inconsistent number of gene extraction before and after')
    return s1qdf,s2qdf


def store_Combed(comdf):
    tree_Chr = {}
    for index,row in comdf.iterrows():
        if row['score'] < 0:
            AB = 'B'
        elif row['score'] > 0:
            AB = 'A'
        else:
            AB = 'NaN'

        if row['seqnames'] not in tree_Chr:
            tree_Chr[row['seqnames']] = IntervalTree()
            tree_Chr[row['seqnames']].addi(row['start'], row['end'], AB)
        else:
            tree_Chr[row['seqnames']].addi(row['start'], row['end'], AB)
    return tree_Chr


def store_genebed(genedf):
    tree_gene = {}
    for index,row in genedf.iterrows():
        if row['seqid'] not in tree_gene:
            tree_gene[row['seqid']] = IntervalTree()
            tree_gene[row['seqid']].addi(row['start'],row['end'],row['ID'])
        else:
            tree_gene[row['seqid']].addi(row['start'], row['end'], row['ID'])
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
    from rich.console import Console
    from rich.traceback import install
    console = Console()
    install()

    s1com = args.s1com
    s2com = args.s2com
    anchors = args.anchors
    s1jcvibed = args.s1bed
    s2jcvibed = args.s2bed
    s1output = args.s1out
    s2output = args.s2out

    s1comdf = read_compartmentfile(s1com)
    s2comdf = read_compartmentfile(s2com)
    s1ancID,s2ancID = read_anchors(anchors)
    s1allgenebed = read_allgenebed(s1jcvibed)
    s2allgenebed = read_allgenebed(s2jcvibed)
    s1qbed,s2qbed = Find_query_genebed(s1ancID,s2ancID,s1allgenebed,s2allgenebed)
    s1tree_gene = store_genebed(s1qbed)
    s2tree_gene = store_genebed(s2qbed)
    s1tree_Chr = store_Combed(s1comdf)
    s2tree_Chr = store_Combed(s2comdf)
    s1ced_t_gene = classifyGene2AB(s1tree_Chr, s1tree_gene)
    s2ced_t_gene = classifyGene2AB(s2tree_Chr, s2tree_gene)
    Output(s1ced_t_gene, s1output)
    Output(s2ced_t_gene, s2output)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Based on the compartment results with the gene ID file on the chromosome, assign each gene to the A/B compartment and count its percentage
        根据compartment结果与染色体上的基因ID文件，给每个基因分配到A/B compartment并统计其占比''',
        usage="python {} -a -o1 ".format(sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,
                                                                            __mail__, __date__, __version__))
    parser.add_argument('-a', '--anchors', required=True, help='Input query anchors file')

    parser.add_argument('-1c', '--s1com', required=True, help='Input species1 compartment result(.csv)')
    parser.add_argument('-2c', '--s2com', required=True, help='Input species2 compartment result(.csv)')

    parser.add_argument('-1b', '--s1bed', required=True, help='Input species1 jcvi gene bed')
    parser.add_argument('-2b', '--s2bed', required=True, help='Input species2 jcvi gene bed')

    parser.add_argument('-1o', '--s1out', required=True, help='Output species1 gene-AB assigned results(.txt)')
    parser.add_argument('-2o', '--s2out', required=True, help='Output species2 gene-AB assigned results(.txt)')
    args = parser.parse_args()
    main(args)
