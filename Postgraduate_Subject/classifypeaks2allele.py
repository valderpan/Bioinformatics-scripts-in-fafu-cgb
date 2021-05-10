#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/5/5


import pandas as pd
from collections import Counter
import math
import numpy as np

def classifyAllelegene(genome_gff3):
    allele_A = [];allele_B = [];allele_C = [];allele_D = [];allgene = []
    with open(genome_gff3,'r') as f:
        lines = (line.strip() for line in f)
        for line in lines:
            if line.startswith('#'):
                continue
            else:
                line_list = line.split('\t')
                if line_list[2] == 'gene':
                    # gene2chr = {}
                    geneid = line_list[8].split(';')[0].split('=')[1]
                    allgene.append(geneid)
                    if line_list[0].endswith('A'):
                        allele_A.append({geneid:line_list[0]})
                    elif line_list[0].endswith('B'):
                        allele_B.append({geneid:line_list[0]})
                    elif line_list[0].endswith('C'):
                        allele_C.append({geneid:line_list[0]})
                    elif line_list[0].endswith('D'):
                        allele_D.append({geneid: line_list[0]})
    print(len(allgene))
    return allgene,allele_A,allele_B,allele_C,allele_D


def read_peak_file(peak_file):
    if peak_file.endswith('.csv'):
        df = pd.read_csv(peak_file)
    elif peak_file.endswith('.xlsx'):
        df = pd.read_excel(peak_file)
    else:
        df = pd.read_table(peak_file,sep='\t')
    return df


def classifypeaks2Chr_allele(peak_df,output_file=0):
    chr_list = peak_df['seqnames'].unique().tolist()
    chr2peaksnum = {}
    chr2peakslen = {}
    for chr in chr_list:
        query_df = peak_df[peak_df['seqnames'] == chr]
        # print(query_df)
        chr2peaksnum[chr] = query_df.shape[0]
        chr2peakslen[chr] = query_df['width'].sum()
    chr_peaksnum = pd.DataFrame([chr2peaksnum]).T
    chr_peakslen = pd.DataFrame([chr2peakslen]).T
    # print({'chr':result})
    # print(chr_peaksnum)
    # print(chr_peakslen)
    result = pd.merge(chr_peaksnum,chr_peakslen,on=chr_peaksnum.index)
    # print(result)
    if output_file:
        result.to_csv('peaks_onchr_distribution.tab',sep='\t',header=False,index=False)

# def classpeaks2allele(peak_file):
#     df = pd.read_excel(peak_file)
#     df_geneid = df['geneId'].tolist()
#     df_geneid_count = dict(Counter(df_geneid))
#     print(df)
#     df_chr = df['seqnames'].unique().tolist()
#     print(df_chr)
#     # print(df_geneid)
#     # print(df_geneid_count)
#     # print(len(df_geneid))
#     # print(len(df_geneid_count))
#     gene_ACR_1 = 0; gene_ACR_2 = 0;gene_ACR_3 = 0;gene_ACR_4 = 0
#     gene_acr_freq_list = []
#     for geneid in df_geneid_count.keys():
#         if df_geneid_count[geneid] not in gene_acr_freq_list:
#             gene_acr_freq_list.append(df_geneid_count[geneid])
#         else:
#             continue
#     print(gene_acr_freq_list)
#     for gene in df_geneid_count.keys():
#         print(gene)
#         break
#
#     # print(df_geneid_count.items())
#     # for num in gene_acr_freq_list:
#         # if df_geneid_count.items()[1] == num:
#         #     gene_ACR_1 += 1
#         # elif df_geneid_count[geneid] == 2:
#         #     gene_ACR_2 += 1
#         # elif df_geneid_count[geneid] == 3:
#         #     gene_ACR_3 += 1
#         # elif df_geneid_count[geneid] == 4:
#         #     gene_ACR_4 += 1
#     gene_acr_frqe = {}


def cal_ACR_length(peak_df):
    total_length = peak_df['width'].sum()
    print('ACR total lenth : {}bp'.format(total_length))
    return total_length


def cal_ACR_ratio(acrlength,genome_size):
    acr_ratio = round(acrlength/genome_size*100,2)
    print('Percentage of total genome :{}%'.format(acr_ratio))


def classifyACR3category(peak_df,output_file = 0,tooutput=0):
    '''
    将所有的ACR分为gACR、dACR、pACR三类
    :param peak_df:
    :return:
    '''
    ACR3cate = []
    for index,row in peak_df.iterrows():
        if row['annotation'] == 'Distal Intergenic':
            ACR3cate.append('dACR')
        elif row['annotation'].startswith('Intron'):
            ACR3cate.append('gACR')
        elif row['annotation'].startswith('Exon'):
            ACR3cate.append('gACR')
        elif row['annotation'] == 'Promoter (<=1kb)':
            ACR3cate.append('pACR')
        elif row['annotation'] == 'Promoter (1-2kb)':
            ACR3cate.append('pACR')
        elif row['annotation'] == 'Promoter (2-3kb)':
            ACR3cate.append('dACR')
        elif row['annotation'] == 'Downstream (<1kb)':
            ACR3cate.append('pACR')
        elif row['annotation'] == 'Downstream (1-2kb)':
            ACR3cate.append('pACR')
        elif row['annotation'] == 'Downstream (2-3kb)':
            ACR3cate.append('dACR')
        elif row['annotation'] == "3' UTR":
            ACR3cate.append("gACR")
        elif row['annotation'] == "5' UTR":
            ACR3cate.append("gACR")
    peak_df['category'] = ACR3cate
    mutate_df = peak_df
    # print(mutate_df)
    if output_file:
        mutate_df.to_excel('SES208_peak_category.xlsx',header=True,index=False)

    #--------所有染色体上3种类别的统计-------
    cate2num = {}
    for cate in mutate_df['category'].unique().tolist():
        query_cate_num = mutate_df[mutate_df['category'] == cate].shape[0]
        cate2num[cate] = query_cate_num
    cate2num_df = pd.DataFrame([cate2num])
    cate2num_df['total'] = cate2num_df.T[0].sum()
    cate2num_df = cate2num_df.T
    # print(cate2num_df)

    #--------每条染色体上3中类别的统计--------
    cate2Chrnum = {}
    for chrid in mutate_df['seqnames'].unique().tolist():
        for cate in mutate_df['category'].unique().tolist():
            query_cate_chr = mutate_df[mutate_df['seqnames'] == chrid]
            query_cate_chr = query_cate_chr[query_cate_chr['category'] == cate]
            query_cate_chrnum = query_cate_chr.shape[0]
            cate2Chrnum[(chrid,cate)] = query_cate_chrnum
    # print(cate2Chrnum)
    if tooutput:
        for key in cate2Chrnum.keys():
            print(key[0],key[1],cate2Chrnum[key],sep='\t')
    # result = pd.DataFrame([cate2Chrnum]).T
    # print(result)
    return  mutate_df


def gene3ACRcategory(peak_mutate_df,allgene,rnaseq_file,tooutput=0):
    category2geneid = {}
    for cate in peak_mutate_df['category'].unique().tolist():
        query_cate_df = peak_mutate_df[peak_mutate_df['category'] == cate]
        category2geneid[cate] = query_cate_df['geneId'].tolist()
    # print(category2geneid.keys())
    category2uniqgene = {}
    for cate in category2geneid.keys():
        # print(cate)
        # print(len(category2geneid[cate]))
        if cate == 'dACR':
            category2uniqgene['only{}'.format(cate)] = [gene for gene in category2geneid[cate] if gene not in category2geneid['gACR'] and gene not in category2geneid['pACR']]
            category2uniqgene['gdACR'] = [gene for gene in category2geneid[cate] if gene in category2geneid['gACR'] and gene not in category2geneid['pACR']]
            category2uniqgene['dpACR'] = [gene for gene in category2geneid[cate] if gene not in category2geneid['gACR'] and gene in category2geneid['pACR']]
            category2uniqgene['gdpACR'] = [gene for gene in category2geneid[cate] if gene in category2geneid['gACR'] and gene in category2geneid['pACR']]
        elif cate == 'gACR':
            category2uniqgene['only{}'.format(cate)] = [gene for gene in category2geneid[cate] if gene not in category2geneid['dACR'] and gene not in category2geneid['pACR']]
            category2uniqgene['gpACR'] = [gene for gene in category2geneid[cate] if gene not in category2geneid['dACR'] and gene in category2geneid['pACR']]
        elif cate == 'pACR':
            category2uniqgene['only{}'.format(cate)] = [gene for gene in category2geneid[cate] if gene not in category2geneid['gACR'] and gene not in category2geneid['dACR']]

    # for i in category2uniqgene.keys():
    #     print(i)
    #     print(len(category2uniqgene[i]))


    # print(len(allgene))
    category2uniqgene['allgene'] = allgene
    # print(len(peak_mutate_df['geneId'].unique().tolist()))
    non_acr = []
    # print()
    for gene in allgene:
        if gene not in peak_mutate_df['geneId'].unique().tolist():
            non_acr.append(gene)
    # print(len(non_acr))
    category2uniqgene['non_ACR'] = non_acr


    rnaseq_df = pd.read_excel(rnaseq_file)
    if tooutput:
        ACR7category_writer = pd.ExcelWriter('ACR7category.xlsx')
        for key in category2uniqgene:
            query_df = rnaseq_df[rnaseq_df['gene_id'].isin(category2uniqgene[key])]
            query_df = query_df[query_df['FPKM'] >=1]
            query_df['FPKM'] = query_df['FPKM'].apply(np.log10)
            query_df.to_excel(ACR7category_writer,sheet_name=key,header=True,index=False)
        ACR7category_writer.close()


def stat_perGeneACRs(peak_df,rnaseq_file,tooutput=0):
    '''
    输出含有不同数目ACR的基因，结合RNA-Seq分析表达情况
    :param peak_df:
    :return:
    '''
    print(peak_df)
    df = peak_df[peak_df['category'] != 'dACR']
    print(df)
    df = df.groupby('geneId')['newanno'].value_counts().unstack()
    # print(df)
    df['ACRs'] = df.apply(lambda x:x.sum(),axis=1)
    df = df.reset_index().iloc[:,[0,-1]]
    df['ACRs'] = df['ACRs'].map(int)
    print(df)
    num2geneid1 = {}
    for num in df['ACRs'].unique().tolist():
        query_df = df[df['ACRs'] == num]
        # print(query_df)
        num2geneid1[num] = query_df['geneId'].unique().tolist()


    up3gene = []
    num2geneid2 = {}
    for i in num2geneid1.keys():
        if i <= 2:
            # print(i)
            # print(len(num2geneid1[i]))
            num2geneid2[i] = num2geneid1[i]
        if i >2:
            up3gene.extend(num2geneid1[i])
    num2geneid2['up3'] = up3gene

    rnaseq_df = pd.read_excel(rnaseq_file)
    if tooutput:
        pergeneACRS_writer = pd.ExcelWriter('pergeneACRs_onlygpACR.xlsx')
        for key in num2geneid2.keys():
            query_num_df = rnaseq_df[rnaseq_df['gene_id'].isin(num2geneid2[key])]
            query_num_df = query_num_df[query_num_df['FPKM'] >= 1]
            query_num_df['FPKM'] = query_num_df['FPKM'].apply(np.log10)
            query_num_df.to_excel(pergeneACRS_writer, sheet_name=str(key), header=True, index=False)
        pergeneACRS_writer.close()


if __name__ == '__main__':
    df = read_peak_file(r'E:\潘浩然\调控元件生信任务\analysis\Sugarcane\peaks\SES208_peak_mutated.xlsx')
    print(df)
    allgene,alleleA,alleleB,alleleC,allelD = classifyAllelegene(r'E:\潘浩然\调控元件生信任务\analysis\callpeak\Sspon.v20190103.gff3')
    classifypeaks2Chr_allele(df)
    acr_length = cal_ACR_length(df)
    cal_ACR_ratio(acr_length,3140615268)
    peak_mutate_df=classifyACR3category(df)
    print(peak_mutate_df)
    # gene3ACRcategory(peak_mutate_df,allgene,r'E:\潘浩然\调控元件生信任务\analysis\Sugarcane\peaks\ses208_mature_leaf_RNAseq.xlsx')
    stat_perGeneACRs(peak_mutate_df,r'E:\潘浩然\调控元件生信任务\analysis\Sugarcane\peaks\ses208_mature_leaf_RNAseq.xlsx',1)