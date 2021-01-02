#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/1/1


import sys
import Fontcolor as Font


def Gain_gmapgene2(gmap_gff3):
    '''
    将gmap.gff3文件存入字典 →→→ {geneID:[tig_list]}
    :param gmap_gff3:
    :return:
    '''
    gene2ctg_dict = {}
    with open(gmap_gff3) as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                line_list = line.strip().split('\t')
                if line_list[2] == 'gene':
                    #Name=Sobic.001G000100→geneID
                    geneID = '.'.join(line_list[8].split(';')[1].split('=')[1].split('.')[0:2])
                    if geneID not in gene2ctg_dict.keys():
                        gene2ctg_dict[geneID] = [line_list[0]]
                    else:
                        gene2ctg_dict[geneID].append(line_list[0])
    return gene2ctg_dict


def parse_refgff3(ref_gff):
    '''
    将refrence_gff3文件存入字典 →→→ {geneID:chrID@start}
    :param ref_gff:
    :return:
    '''
    gene2start_dict = {}
    with open(ref_gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                line_list = line.strip().split()
                if line_list[2] == 'gene':
                    #Name=Sobic.001G000100→geneID
                    geneID = line_list[8].split(';')[1].split('=')[1]
                    gene2start_dict[geneID] = line_list[0]+'@'+line_list[3]
    return gene2start_dict


def match_tig2gene(gene2ctgDict,gene2startDict):
    '''
    匹配并输出
    :param gene2ctgDict:
    :param gene2startDict:
    :return:
    '''
    for gene in sorted(gene2ctgDict.keys()):
        print('\t'.join(gene2startDict[gene].split('@'))+'\t'+'\t'.join(set(gene2ctgdict[gene])))


if __name__ == '__main__':
    if len(sys.argv) != 3:
        Font.Scripts_tip('Parse the results of gmap to generate Allele.ctg.Table')
        print('Usage:')
        print('\tpython {0} <gmap.gff3> <referenece.gff3> > <Allele.ctg.Table>'.format(sys.argv[0]))
    else:
        gene2ctgdict = Gain_gmapgene2(sys.argv[1])
        gene2startdict = parse_refgff3(sys.argv[1])
        match_tig2gene(gene2ctgdict,gene2startdict)
