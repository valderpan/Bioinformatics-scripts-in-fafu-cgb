#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2020/11/6

'''
此脚本用来对等位表进行处理，分别将含有1/2/3/4个等位的基因ID过滤后输出
'''

import pandas as pd
import argparse
import sys

__author__ = 'Haoran Pan'
__mail__ = 'haoran_pan@qq.com'
__date__ = '20201106'
__version__ = 'v1.1'



def allele_selection(alleletable):
    '''
    根据等位表将含有1/2/3/4个等位的基因进行筛选并输出
    :param alleletable: genome alleleTable
    :return: allele1,allele2,allele3,allele4
    '''
    alleleTable = pd.read_csv(alleletable,skiprows=12)
    alleleTable = alleleTable.iloc[:,[2,7,12,17]]
    alleleTable.rename(columns={'Representive Gene model':'Representive_Gene_model','Gene ID':'alleleA_ID',
                                'Gene ID.1':'alleleB_ID','Gene ID.2':'alleleC_ID','Gene ID.3':'alleleD_ID'},inplace=True)
    # print(alleleTable)

    allele_1concat = [];allele_2concat = [];allele_3concat = [];allele_4concat = []
    for index,row in alleleTable.iterrows():
        row_count = list(alleleTable.loc[index].value_counts())
        if len(row_count) == 1:
            allele_1concat.append(alleleTable.loc[[index], :])
        elif len(row_count) == 2:
            allele_2concat.append(alleleTable.loc[[index], :])
        elif len(row_count) == 3:
            allele_3concat.append(alleleTable.loc[[index], :])
        elif len(row_count) == 4:
            allele_4concat.append(alleleTable.loc[[index], :])
    allele_1 = pd.concat(allele_1concat)
    allele_2 = pd.concat(allele_2concat)
    allele_3 = pd.concat(allele_3concat)
    allele_4 = pd.concat(allele_4concat)
    return allele_1,allele_2,allele_3,allele_4


def clean_diff_allele(allele_1,allele_2,allele_3,allele_4):
    '''
    将筛选后的含有1/2/3/4等位的等位基因行进行清洗，去掉‘|’后边的冗余等位并输出
    :param allele_1: Allele table containing 1 allele
    :param allele_2: Allele table containing 2 allele
    :param allele_3: Allele table containing 3 allele
    :param allele_4: Allele table containing 4 allele
    :return:Allele table after filter
    '''
    allele_1_filter = [];allele_2_filter = [];allele_3_filter = [];allele_4_filter = []
    for index,row in allele_1.iterrows():
        row_list = row.tolist()
        row_concat = []
        for i in row_list:
            if '|' in str(i):
                row_concat.append(str(i).split('|')[0])
            else:
                row_concat.append(i)
        # print(row_concat)
        allele_1_filter.append(row_concat)
    allele1 = pd.DataFrame(allele_1_filter)
    # print(allele1)
    for index,row in allele_2.iterrows():
        row_list = row.tolist()
        row_concat = []
        for i in row_list:
            if '|' in str(i):
                row_concat.append(str(i).split('|')[0])
            else:
                row_concat.append(i)
        # print(row_concat)
        allele_2_filter.append(row_concat)
    allele2 = pd.DataFrame(allele_2_filter)

    for index,row in allele_3.iterrows():
        row_list = row.tolist()
        row_concat = []
        for i in row_list:
            if '|' in str(i):
                row_concat.append(str(i).split('|')[0])
            else:
                row_concat.append(i)
        # print(row_concat)
        allele_3_filter.append(row_concat)
    allele3 = pd.DataFrame(allele_3_filter)

    for index,row in allele_4.iterrows():
        row_list = row.tolist()
        row_concat = []
        for i in row_list:
            if '|' in str(i):
                row_concat.append(str(i).split('|')[0])
            else:
                row_concat.append(i)
        # print(row_concat)
        allele_4_filter.append(row_concat)
    allele4 = pd.DataFrame(allele_4_filter)

    return allele1,allele2,allele3,allele4


def gain_AlleleRepresentiveGene(allele1,allele2,allele3,allele4):
    '''
    将不同数目等位基因的具有代表性的(最前面)的基因输出出来，输出geneIDlist
    :param allele1:
    :param allele2:
    :param allele3:
    :param allele4:
    :return:
    '''
    allele1_Repgene = [];allele2_Repgene = [];allele3_Repgene = [];allele4_Repgene = []
    for index,row in allele1.iterrows():
        if type(row[0]) != float:
            allele1_Repgene.append(row[0])
        else:
            if type(row[1]) != float:
                allele1_Repgene.append(row[1])
            else:
                if type(row[2]) != float:
                    allele1_Repgene.append(row[2])
                else:
                    if type(row[3]) != float:
                        allele1_Repgene.append(row[3])
    for index,row in allele2.iterrows():
        if type(row[0]) != float:
            allele2_Repgene.append(row[0])
        else:
            if type(row[1]) != float:
                allele2_Repgene.append(row[1])
            else:
                if type(row[2]) != float:
                    allele2_Repgene.append(row[2])
                else:
                    if type(row[3]) != float:
                        allele2_Repgene.append(row[3])
    for index,row in allele3.iterrows():
        if type(row[0]) != float:
            allele3_Repgene.append(row[0])
        else:
            if type(row[1]) != float:
                allele3_Repgene.append(row[1])
            else:
                if type(row[2]) != float:
                    allele3_Repgene.append(row[2])
                else:
                    if type(row[3]) != float:
                        allele3_Repgene.append(row[3])
    for index,row in allele4.iterrows():
        if type(row[0]) != float:
            allele4_Repgene.append(row[0])
        else:
            if type(row[1]) != float:
                allele4_Repgene.append(row[1])
            else:
                if type(row[2]) != float:
                    allele4_Repgene.append(row[2])
                else:
                    if type(row[3]) != float:
                        allele4_Repgene.append(row[3])
    return allele1_Repgene,allele2_Repgene,allele3_Repgene,allele4_Repgene


def gain_AlleleEveryGene(allele1,allele2,allele3,allele4):
    '''
    得到1/2/3/4等位表中所有的geneID
    :param allele1:Allele table containing 1 allele after filter
    :param allele2:Allele table containing 2 allele after filter
    :param allele3:Allele table containing 3 allele after filter
    :param allele4:Allele table containing 4 allele after filter
    :return:Allele geneID
    '''
    allele1_allgene = [];allele2_allgene = [];allele3_allgene = [];allele4_allgene = []

    for j in range(0,4):
        for i in allele1.iloc[:,j].tolist():
            if type(i) != float :
                allele1_allgene.append(i)

    for j in range(0,4):
        for i in allele2.iloc[:,j].tolist():
            if type(i) != float :
                allele2_allgene.append(i)

    for j in range(0,4):
        for i in allele3.iloc[:,j].tolist():
            if type(i) != float :
                allele3_allgene.append(i)

    for j in range(0,4):
        for i in allele4.iloc[:,j].tolist():
            if type(i) != float :
                allele4_allgene.append(i)

    return allele1_allgene,allele2_allgene,allele3_allgene,allele4_allgene


def output_AlleleGene(allele1gene,allele2gene,allele3gene,allele4gene,allele1out,allele2out,allele3out,allele4out):
    '''
    输出1/2/3/4等位表中所有的geneID
    :param allele1gene: Allele geneID in allele table which contains 1 allele
    :param allele2gene: Allele geneID in allele table which contains 2 allele
    :param allele3gene: Allele geneID in allele table which contains 3 allele
    :param allele4gene: Allele geneID in allele table which contains 4 allele
    :param allele1out: Output allele1gene
    :param allele2out: Output allele2gene
    :param allele3out: Output allele3gene
    :param allele4out: Output allele4gene
    :return: file output
    '''
    with open(allele1out,'w') as f1:
        for i in allele1gene:
            f1.write(str(i)+'\n')
    with open(allele2out,'w') as f2:
        for i in allele2gene:
            f2.write(str(i)+'\n')
    with open(allele3out,'w') as f3:
        for i in allele3gene:
            f3.write(str(i)+'\n')
    with open(allele4out,'w') as f4:
        for i in allele4gene:
            f4.write(str(i)+'\n')


def main(args):
    if args.alleletable and args.mode == 'all':
        genome_alleletable = args.alleletable
        allele1out = args.allele1output
        allele2out = args.allele2output
        allele3out = args.allele3output
        allele4out = args.allele4output
        allele_1, allele_2, allele_3, allele_4 = allele_selection(genome_alleletable)
        allele1, allele2, allele3, allele4 = clean_diff_allele(allele_1, allele_2, allele_3, allele_4)
        allele1_allGene,allele2_allGene,allele3_allGene,allele4_allGene=gain_AlleleEveryGene(allele1, allele2, allele3, allele4)
        output_AlleleGene(allele1_allGene,allele2_allGene,allele3_allGene,allele4_allGene,allele1out,allele2out,allele3out,allele4out)
    elif args.alleletable and args.mode == 'representive':
        genome_alleletable = args.alleletable
        allele1out = args.allele1output
        allele2out = args.allele2output
        allele3out = args.allele3output
        allele4out = args.allele4output
        allele_1, allele_2, allele_3, allele_4 = allele_selection(genome_alleletable)
        allele1, allele2, allele3, allele4 = clean_diff_allele(allele_1, allele_2, allele_3, allele_4)
        allele1_RepGene,allele2_RepGene,allele3_RepGene,allele4_RepGene=gain_AlleleRepresentiveGene(allele1,allele2,allele3,allele4)
        output_AlleleGene(allele1_RepGene, allele2_RepGene, allele3_RepGene, allele4_RepGene, allele1out, allele2out,allele3out, allele4out)
    else:
        print('Please run "python {} -h"'.format(sys.argv[0]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''gain allele geneID''',
        conflict_handler='resolve',
        usage="python {} -a alleletable.csv --mode -o1 allele_1.txt -o2 allele_2.txt -o3 allele_3.txt -o4 allele_4.txt".format(sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__, __mail__, __date__,
                                                                            __version__))
    parser.add_argument('-a', '--alleletable', required=True, help='Input Genome Assembly Allele Table(.csv)')
    parser.add_argument('--mode', required=True, type=str,help='Please choose to output all genes or representative genes [all or representive]')
    parser.add_argument('-o1', '--allele1output', required=True, help='Output allele geneID contains 1 allele')
    parser.add_argument('-o2', '--allele2output', required=True, help='Output allele geneID contains 2 allele')
    parser.add_argument('-o3', '--allele3output', required=True, help='Output allele geneID contains 3 allele')
    parser.add_argument('-o4', '--allele4output', required=True, help='Output allele geneID contains 4 allele')
    args = parser.parse_args()
    main(args)