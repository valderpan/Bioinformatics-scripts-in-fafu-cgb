# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-


import os
import subprocess
from datetime import datetime
import re
import pandas as pd
import argparse

#Description : This script is used to convert RSEM quantitative results into a gene expression matrix
#Usage : python rsem2gene_exp_df.py path1 path2
#Note : path1 is the parent directory where rsem output files are stored ; path2 is the parent directory of fpkm-1 and fpkm-2

Folder_name_list = []

def run_command(cmd):
    print(cmd)
    return_code = subprocess.call(cmd,shell=True)
    if return_code != 0 :
        print('ERROR:[{2}] Return code {0} when running the following command : {1} '.format(return_code, cmd,datetime.now()))

def Get_folder_name(folder_path):
    path1 = folder_path
    path_folder = os.listdir(path1)
    for folder in path_folder:
        Folder_name_list.append(folder)

    os.mkdir(path2 + 'fpkm-1')
    os.mkdir(path2 + 'fpkm-2')
    return Folder_name_list

# Step 1: Rename and copy Rsem.genes.result under each sample-ID directory to fpkm-1 directory
def rename_and_copy(folder_path):
    for folder in Folder_name_list:
        os.chdir(folder_path + folder)
        cmd1 = 'mv ' + 'RSEM.genes.results ' + folder+'.genes.results'
        run_command(cmd1)
        cmd2 = 'cp ' + folder +'.genes.results' + ' ../../fpkm-1/'
        run_command(cmd2)


# Step 2: Extract the first and seventh columns from each file in the fpkm-1 directory and enter them into the fpkm-2 directory
def gain_IDandFPKM(path):
    fpkm_1_path = path  + 'fpkm-1'
    fpkm_1_list = os.listdir(fpkm_1_path)
    os.chdir(fpkm_1_path)
    for i in fpkm_1_list:
        cmd3 = 'cat ' + i + '|' + 'cut -f 1,7' + ' >' + i+'.fpkm.result'
        run_command(cmd3)
        cmd4 = 'mv ' + i + '.fpkm.result ' + '../fpkm-2/'
        run_command(cmd4)


# Step 3: Merge each file in the fpkm-2 directory into an expression matrix and output the matrix
def merge_files_to_matrix(path,keyword):
    fpkm_2_path = path + 'fpkm-2'
    os.chdir(fpkm_2_path)
    fpkm_file_list = os.listdir('../')
    fpkm_file = []
    for i in fpkm_file_list:
        # j = re.findall('([a-zA-Z0-9\-]+)\.',i)
        j = re.findall(keyword,i)
        fpkm_file.append(j[0])


    first_data = pd.read_csv(fpkm_file_list[0],sep = '\t',names=['gene_id',fpkm_file[0]])
    first_data = first_data.drop(first_data.index[0])
    file_dict = {}
    # print(first_data)

    for i in range(1,len(fpkm_file_list)):
        file_dict[i] = pd.read_csv(fpkm_file_list[i],sep = '\t',names=['gene_id',fpkm_file[i]])
        if i < len(fpkm_file_list):
            first_data = pd.merge(first_data,file_dict[i],on='gene_id')

    first_data.to_excel('fpkm_dataframe.xlsx',header=True,index=False)


# Step 4: Output the matrix file and delete the generated intermediate directory
def Delete_the_created_files(path,folder_path):
    cmd5 = 'mv ' + 'fpkm_dataframe.xlsx' + ' ../'
    run_command(cmd5)
    os.chdir(path)
    cmd6 = 'rm -rf ' + 'fpkm-1'
    run_command(cmd6)
    cmd7 = 'rm -rf ' + 'fpkm-2'
    run_command(cmd7)
    for folder in Folder_name_list:
        os.chdir(folder_path + folder)
        cmd8 = 'mv ' + folder+'.genes.results '+ 'RSEM.genes.results'
        run_command(cmd8)

def main(args):
    folder_path = args.path1
    global path2
    path2 = args.path2
    keyword = args.key

    Get_folder_name(folder_path)
    print('Step 1: Rename and copy Rsem.genes.result under each sample-ID directory to fpkm-1 directory.\n')
    rename_and_copy(folder_path)

    print('Step 2: Extract the first and seventh columns from each file in the fpkm-1 directory and enter them into the fpkm-2 director.\n')
    gain_IDandFPKM(path2)

    print('Step 3: Merge each file in the fpkm-2 directory into an expression matrix and output the matrix.\n')
    merge_files_to_matrix(path2,keyword)

    print('Step 4: Output the matrix file and delete the generated intermediate directory.\n')
    Delete_the_created_files(path2,folder_path)

    print('All tasks were finished!\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='rsem2gene_exp_df.py',
        formatter_class=argparse.RawTextHelpFormatter,
        description='''This script is used to convert RSEM quantitative results into a gene expression matrix''')
    parser.add_argument('-path1',required=True,help='path1 is the parent directory where rsem output files are stored')
    parser.add_argument('-path2',required=True,help='path2 is the parent directory of fpkm-1 and fpkm-2')
    parser.add_argument('-key',required=True,help='Please enter regular expression keywords')
    args = parser.parse_args()
    main(args)