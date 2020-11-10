#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2020/11/6

import pandas as pd
from plotnine import *
import sys
def get_geneCRE_base_geneID(geneID,allgenomeCRE,output):
    allgeneCRE = pd.read_table(allgenomeCRE,sep='\t')
    allgeneCRE['GeneID'] = allgeneCRE['GeneID'].str[:-1]
    gene_df = pd.read_table(geneID,header=None)
    gene_list = gene_df.iloc[:,0].tolist()
    query_genedf = allgeneCRE[allgeneCRE['GeneID'].isin(gene_list)]
    query_genedf = query_genedf.iloc[:, [0, 1, -1]]
    query_genedf['Count'] = query_genedf['Position'].str.count(',') + 1
    query_genedf = query_genedf.groupby(['GeneID','Site Name'])['Count'].sum().reset_index()
    query_genedf = query_genedf.groupby(by='Site Name').agg({"Count":sum}).reset_index()
    query_genedf = query_genedf.sort_values(by='Count',ascending=False)
    # query_genedf.to_csv(output,sep='\t',header=True,index=False)
    top20 = query_genedf.iloc[0:20]
    top20.to_csv(output, sep='\t', header=True, index=False)
# def plot_pie_for_geneCRE(query_gene_df):
#     tmp = (
#         ggplot(query_gene_df,aes(x='',y='Count',fill='Site Name'))+
#         geom_bar(stat='identity',width=1,color='white')+
#         coolrd('y',start=0)
#     )
#     tmp.save('test.pdf',width=10,height=15)

if __name__ == '__main__':
    # get_geneCRE_base_geneID(r'E:\潘浩然\调控元件生信任务\数据\allele_CRE_compare\allele1.txt',
    #                         r'E:\潘浩然\调控元件生信任务\数据\Sspon_all_genome_CREv1101.tab',
    #                         r'E:\潘浩然\调控元件生信任务\数据\allele_CRE_compare\allele1CRE_infotop20.tab')
    # plot_pie_for_geneCRE(query_gene_df)
    get_geneCRE_base_geneID(sys.argv[1],sys.argv[2],sys.argv[3])