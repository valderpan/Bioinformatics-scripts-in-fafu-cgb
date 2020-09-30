#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# @Time : 2020/8/20 16:37


import pandas as pd
import os
import argparse

def concat_eggNoG_results(eggNoG_resultspath,output_file):
    os.chdir(eggNoG_resultspath)
    file_list = os.listdir(eggNoG_resultspath)
    concat_df_list = []
    for file in file_list:
        df = pd.read_csv(file,sep = '\t',comment='#',header=None)
        # print(df)
        concat_df_list.append(df)
    result = pd.concat(concat_df_list)
    print(result)
    result.to_csv(output_file,sep='\t',header=False,index=False)
# concat_eggNoG_results(r'C:\Users\dell\Desktop\脚本测试\20200819 clusterProfiler + eggNOG\Nepal_all_genome')
#                     r'C:\Users\dell\Desktop\脚本测试\20200819 clusterProfiler + eggNOG\test.tab')



def match_GO_to_Description(go_tb,Go_anno,output_file):
    go_database = pd.read_table(go_tb,sep='\t')
    go_des_dict = go_database.set_index('GO').to_dict()['Description']
    GOannotation = pd.read_table(Go_anno,sep='\t')
    GOannotation['description'] = GOannotation['GO'].map(go_des_dict)
    print(GOannotation)
    GOannotation.to_csv(output_file,sep='\t',header=True,index=False)


# match_GO_to_Description(r'C:\Users\dell\Desktop\脚本测试\20200819 clusterProfiler + eggNOG\go.tb',
#                         r'C:\Users\dell\Desktop\脚本测试\20200819 clusterProfiler + eggNOG\GOannotation.tsv',
#                         r'C:\Users\dell\Desktop\脚本测试\20200819 clusterProfiler + eggNOG\go_anno.tab')
def main(args):
    if args.path and args.result :
        eggNoG_file_path = args.path
        cleaned_result = args.result
        concat_eggNoG_results(eggNoG_file_path,cleaned_result)
    if args.tb and args.anno and args.output:
        go_tb = args.tb
        go_anno = args.anno
        output = args.output
        match_GO_to_Description(go_tb,go_anno,output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='clean_eggNOG_and_match_GO.py',
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Clean the output results of eggNOG and GO annotation files''')
    parser.add_argument('-p','--path',help='Input the path to save the eggNOG file(.annotations)')
    parser.add_argument('-r','--result',help='Output cleaned eggNOG result file(.annotations)')
    parser.add_argument('-t','--tb',help='Input go.tb')
    parser.add_argument('-a','--anno',help='Input GOannotation.tsv')
    parser.add_argument('-o','--output',help='Output the file to clusterprofiler(.tab)')
    args = parser.parse_args()
    main(args)