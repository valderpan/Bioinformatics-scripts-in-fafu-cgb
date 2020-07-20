# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

import pandas as pd
import argparse

def Find_allele_table_with_4_genes_or_more(allele_file,result_file):
    '''
    step1 根据等位表找出含有三个或者四个等位的基因
    example(allele_file:Sspon.alleleTable20180810.xlsx;
            output_file:LA_allele_4ormore.xlsx)
    注意：运行前需要对等位表进行处理，去掉多余的表头及描述列并保存为excel格式！
    '''
    assembly_allele_df = pd.read_excel(allele_file,header=0)
    print(assembly_allele_df)

    concat_row_list = []
    for index,row in assembly_allele_df.iterrows():
        row_count = list(assembly_allele_df.loc[index].value_counts())
        # print(len(row_count))
        if len(row_count)>4:
            concat_row_list.append(assembly_allele_df.loc[[index],:])

    allele_4ormore = pd.concat(concat_row_list)
    print(allele_4ormore)
    allele_4ormore.to_excel(result_file, header=True, index=False)

def Find_min_KS(Ks_file,allele_result_file,output_file):
    '''
    step2 根据上一步得到的符合条件的等位表，将每行等位中最小基因的ks值筛选出来
    example(Ks_file:LA_Os_KaKs_caculate.result KA/KS计算结果即可
            allele_result_file:\LA_allele_4ormore.xlsx 上一步结果
            output_file:LA_OS_allele_minKs.result)
    '''
    Species_ks = pd.read_table(Ks_file,sep='\t')
    Ks_df = Species_ks['name'].str.split(';',expand=True)
    Ks_df['dS-yn'] = Species_ks['dS-yn']
    Ks_df['dN-yn'] = Species_ks['dN-yn']
    Ks_df['dS-ng'] = Species_ks['dS-ng']
    Ks_df['dN-ng'] = Species_ks['dN-ng']
    Ks_df.rename(columns={0:'Species1_ID',1:'Species2_ID'},inplace=True)
    Ks_df.index = Ks_df['Species1_ID']
    Ks_df = Ks_df.iloc[:,1:]
    print(Ks_df)
    result_concatrow_list = []
    allele_df = pd.read_excel(allele_result_file,header=0)
    print(allele_df)
    for index,row in allele_df.iterrows():
        row_list = row.to_list()
        concat_allelerow_list = []
        # print(row_list)
        for i in row_list:
            if '|' not in str(i):
                if i in Ks_df.index:
                    # print(Ks_df.loc[[i],:])
                    concat_allelerow_list.append(Ks_df.loc[[i], :])
            else:
                first_geneID = str(i).split('|')[0]
                if first_geneID in Ks_df.index:
                    concat_allelerow_list.append(Ks_df.loc[[first_geneID], :])
        if not len(concat_allelerow_list) == 0:
            allelerow_ks = pd.concat(concat_allelerow_list)
            # print(allelerow_ks)
            allelerow_ks.sort_values(by='dS-ng', ascending=True, inplace=True)
            # print(allelerow_ks)
            result_concatrow_list.append(allelerow_ks.head(1))
    result = pd.concat(result_concatrow_list)
    result = result.reset_index()
    result['name'] = result['Species1_ID']+';'+result['Species1_ID']
    result = result.iloc[:,[6,2,3,4,5]]
    print(result)
    result.to_csv(output_file,sep='\t',header=True,index=False)

def main(args):
    if args.assembly_allele:
        allele_file = args.assembly_allele
        result_file = args.allele_result
        Find_allele_table_with_4_genes_or_more(allele_file,result_file)
    if args.ksfile and args.alleleresult:
        Ks_file = args.ksfile
        allele_result_file = args.alleleresult
        output_file = args.output
        Find_min_KS(Ks_file,allele_result_file,output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='find_minKS_in_allele.py',
        description='''Find the gene with the Minimum KS value between alleles''')
    parser.add_argument('-aa','--assembly_allele',help='Input Allele table after genome assembly()(.xlsx)')
    parser.add_argument('-ar','--allele_result',help='Output filtered allele table(.xlsx) ')
    parser.add_argument('-k','--ksfile',help='Input ks file(.ks.result)')
    parser.add_argument('-a','--alleleresult',help='Input the filtered allele table obtained from the previous step(.xlsx)')
    parser.add_argument('-o','--output',help='Output the final result(.ks.result)')
    args = parser.parse_args()
    main(args)
