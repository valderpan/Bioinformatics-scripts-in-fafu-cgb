#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# @Time : 2020/10/7 17:09


'''
需求：找到割手密基因组中含有三组和四组等位的基因,将每组等位中的3和4个基因的CRE输出
'''

import pandas as pd
import argparse
import sys


# pd.set_option('display.max_columns',5)

def Find_allele_table_with_3or4_genes_or_more(allele_file):
    '''
    step1 根据等位表找出含有三个和四个等位的基因
    example(allele_file:Sspon.alleleTable20180810.xlsx;
            output_file:LA_allele_4ormore.xlsx)
    注意：运行前需要对等位表进行处理，去掉多余的表头及描述列并保存为excel格式！
    '''
    assembly_allele_df = pd.read_excel(allele_file,header=0)

    allele_3concat_row_list = [];allele_4concat_row_list = []
    for index,row in assembly_allele_df.iterrows():
        row_count = list(assembly_allele_df.loc[index].value_counts())
        # print(len(row_count))
        if len(row_count)==3:
            allele_3concat_row_list.append(assembly_allele_df.loc[[index],:])
        elif len(row_count) == 4:
            allele_4concat_row_list.append(assembly_allele_df.loc[[index],:])
    allele_3 = pd.concat(allele_3concat_row_list)
    allele_4 = pd.concat(allele_4concat_row_list)
    # allele_3ormore.to_excel(result_file, header=True, index=False)
    return allele_3,allele_4

def get_3or4_allele_gene_id(allele_file,allgenome_CRE_file,output_path):
    '''
    step2:在基因中将‘|’后面的基因删去，只保留第一个geneid;将整理后的矩阵中 每一行的3或4个基因的CRE矩阵输出
    example(allele_file:Sspon.alleleTable20180810.xlsx,
            allgenome_CRE_file:Sspon_all_genome_CRE.tab,
            output_path:output)
    '''
    df3,df4 = Find_allele_table_with_3or4_genes_or_more(allele_file)
    # df = pd.read_excel(allele_file)
    df3concat_list = [];df4concat_list = []
    for index,row in df3.iterrows():
        row_list = row.tolist()
        row_concat = []
        for i in row_list:
            if '|' in str(i):
                row_concat.append(str(i).split('|')[0])
            else:
                row_concat.append(i)
        df3concat_list.append(row_concat)
    result_df3 = pd.DataFrame(df3concat_list)
    # print(result)
    genome_CRE_df = pd.read_table(allgenome_CRE_file,sep='\t')
    # print(genome_CRE_df)

    for index,row in result_df3.iterrows():
        query_df = genome_CRE_df[genome_CRE_df['GeneID'].str[:-1].isin(row.tolist())]
        if len(query_df['GeneID'].unique()) == 3:
            query_df.to_csv(output_path+'/'+'row_'+str(index)+'.allele3_3.tab',sep='\t',header=True,index=False)
        elif len(query_df['GeneID'].unique()) == 2:
            query_df.to_csv(output_path + '/' + 'row_' + str(index) + '.allele3_2.tab', sep='\t', header=True,
                            index=False)
        elif len(query_df['GeneID'].unique()) == 1:
            query_df.to_csv(output_path + '/' + 'row_' + str(index) + '.allele3_1.tab', sep='\t', header=True,
                            index=False)

    for index,row in df4.iterrows():
        row_list = row.tolist()
        row_concat = []
        for i in row_list:
            if '|' in str(i):
                row_concat.append(str(i).split('|')[0])
            else:
                row_concat.append(i)
        df4concat_list.append(row_concat)
    result_df4 = pd.DataFrame(df4concat_list)
    # print(result)
    for index,row in result_df4.iterrows():
        query_df = genome_CRE_df[genome_CRE_df['GeneID'].str[:-1].isin(row.tolist())]
        if len(query_df['GeneID'].unique()) == 4:
            query_df.to_csv(output_path+'/'+'row_'+str(index)+'.allele4_4.tab',sep='\t',header=True,index=False)
        elif len(query_df['GeneID'].unique()) == 3:
            query_df.to_csv(output_path + '/' + 'row_' + str(index) + '.allele4_3.tab', sep='\t', header=True,
                            index=False)
        elif len(query_df['GeneID'].unique()) == 2:
            query_df.to_csv(output_path + '/' + 'row_' + str(index) + '.allele4_2.tab', sep='\t', header=True,
                            index=False)
        elif len(query_df['GeneID'].unique()) == 1:
            query_df.to_csv(output_path + '/' + 'row_' + str(index) + '.allele4_1.tab', sep='\t', header=True,
                            index=False)

def main(args):
    allele_file = args.allele
    allgenome_CRE = args.crematrix
    output_path = args.output
    get_3or4_allele_gene_id(allele_file,allgenome_CRE,output_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Find genes with three or four sets of alleles in the polyploid genome, and output the CREs of 3 or 4 genes in each set of alleles''')
    parser.add_argument('-a', '--allele', required=True, help='Input allele table')
    parser.add_argument('-c', '--crematrix', required=True, help='Import the genome-wide CRE file')
    parser.add_argument('-o', '--output', required=True, help='Output directory path')
    args = parser.parse_args()
    main(args)


