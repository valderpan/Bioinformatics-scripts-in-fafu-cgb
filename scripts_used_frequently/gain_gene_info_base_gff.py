#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# @Time : 2020/9/25 8:32

import re
import pandas as pd
import argparse

def gain_gene_position(gff_file,type,chrid):
    gff_df = pd.read_table(gff_file,comment='#',names=['Chrid','Source','Type','Start','End','Score','Strand','Phase','Attributes'])
    gff_df = gff_df[gff_df['Type'] == type]
    gff_df = gff_df.iloc[:,[0,3,4,8]]
    gff_df['Attributes'] = gff_df['Attributes'].str.split(';',expand=True).iloc[:,0].str.split('=',expand=True).iloc[:,1]
    query_chr_list = list(chrid)
    print(query_chr_list)
    print(query_chr_list[0].split(','))
    if len(query_chr_list[0].split(','))  == 1:
        gff_df = gff_df[gff_df['Chrid'].isin(query_chr_list[0].split(','))]
        chrID =  re.findall('[0-9a-zA-Z]+',gff_df.iloc[1,:]['Chrid'])[0]
        gff_df.to_excel(chrID+'_'+type+'.xlsx',header=False,index=False)
    elif len(query_chr_list[0].split(',')) > 1:
        for i in query_chr_list[0].split(','):
            concat_list = []
            for index,row in gff_df.iterrows():
                if gff_df.loc[index]['Chrid'] == i:
                    concat_list.append(gff_df.loc[[index],:])
            result = pd.concat(concat_list)
            chrID = re.findall('[0-9a-zA-Z]+', result.iloc[1, :]['Chrid'])[0]
            result.to_excel(chrID+'_'+type+'.xlsx',header=False,index=False)

def main(args):
    gff_file = args.gff
    Type = args.type
    Chrid = args.chrid
    gain_gene_position(gff_file,Type,Chrid)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='test.py')
    parser.add_argument('-g', '--gff', required=True, help='Input gff3 file')
    parser.add_argument('-t', '--type', required=True, help='Select the type of data to extract(gene,CDS,exon,five_prime_UTR,mRNA,three_prime_UTR)')
    parser.add_argument('-c', '--chrid', nargs='+',required=True, help='Choose to extract chromosome name')
    args = parser.parse_args()
    main(args)
