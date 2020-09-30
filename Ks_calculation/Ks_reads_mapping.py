# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

import pandas as pd
import reads_mapping as rm
import os
import re
import argparse

file_dict = {}
pd.set_option('display.max_columns',10)
def merge_ks_reads_mapping_df(ks_file_path,Upper_limit,output_file):
    os.chdir(ks_file_path)
    ks_file_list = os.listdir('./')
    # print(ks_file_list)
    ks_file_name_list = []
    for i in ks_file_list:
        j = re.findall('([A-Za-z0-9\_]+)\.',i)
        ks_file_name_list.append(j)
    first_file = rm.reads_mapping_calcultate(ks_file_list[0],ks_file_name_list[0],Upper_limit)
    # print(first_file)
    for i in range(1,len(ks_file_list)):
        file_dict[i] = rm.reads_mapping_calcultate(ks_file_list[i],ks_file_name_list[i],Upper_limit)
        # print(file_dict[i])
        if i < len(ks_file_list):
            first_file = pd.merge(first_file,file_dict[i],on='index',how='outer')
    # print(first_file)
    choose_columns_number =  []
    for i in range(0,first_file.shape[1],2):
        choose_columns_number.append(i)
        # print(i)
    output_ks_file = first_file.iloc[:,choose_columns_number]
    # print(output_ks_file)
    output_ks_file = output_ks_file.fillna(0)
    output_ks_file = output_ks_file.sort_values(by='index')
    # print(output_ks_file)
    output_ks_file.to_excel(output_file,header=True,index=False)


def main(args):
    Ks_file_path = args.path
    Upper_limit = args.limit
    output_file = args.output
    merge_ks_reads_mapping_df(Ks_file_path,Upper_limit,output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='Ks_reads_mapping.py',
        description='''Calculate the percentage of gene pairs with different synonymous substitution ratio''')
    parser.add_argument('-p', '--path', required=True, help='Path to store ks file')
    parser.add_argument('-l', '--limit', required=True, type=float,help='Upper limit of ks reads mapping')
    parser.add_argument('-o', '--output', required=True,help='The name of the output file(.xlsx)')
    args = parser.parse_args()
    main(args)
