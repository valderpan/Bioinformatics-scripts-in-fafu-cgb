# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

'''
RNA-Seq计算表达量
预先创建两个文件夹，一个是01-ref 用来存放参考基因组
'''
import os
import subprocess
import datetime
import re
import argparse

reads_list = []
reads_ID_list = []
def Get_sequencing_data(reads_path):
    all_reads = os.listdir(reads_path)
    for read in all_reads:
        read = read.strip()
        read_ID = re.findall('([a-z0-9A-Z\-\_]+)\_',read)   #Modifiable
        reads_list.append(read_ID)
    for ID in reads_list:
        if ID not in reads_ID_list:
            reads_ID_list.append(ID)
    return reads_ID_list

def run_command(cmd):
    print(cmd)
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
        print("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd,datetime.datetime.now()))

def RSEM_calculate_FPKM(gene_fasta,reads_path):
    os.mkdir('rsem_result')
    reads_ID_str = [''.join(map(str, sub_list)) for sub_list in reads_ID_list]
    for read in reads_ID_str:
        # Modifiable
        cmd = 'align_and_estimate_abundance.pl --transcripts '+gene_fasta+' --seqType fq '+ '--left '+reads_path+read+'_1.clean.fq.gz '+ '--right ' +reads_path+read+'_2.clean.fq.gz '+ '--est_method RSEM '+'--output_dir ' + './rsem_result/'+read + ' --aln_method bowtie2 --thread_count 20 --prep_reference'
        run_command(cmd)

def main(args):
    ref_genome = args.genome
    reads_path = args.path
    Get_sequencing_data(reads_path)
    RSEM_calculate_FPKM(ref_genome,reads_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='RNA-Seq_RSEM.py',
        formatter_class=argparse.RawTextHelpFormatter,
        description='''This script is used to execute RNA-Seq pipeline to used to calculate gene expression values''')
    parser.add_argument(
        '-g', '--genome',
        required=True,
        help='Please specify a reference genome')
    parser.add_argument(
        '-p', '--path',
        required=True,
        help='Please specify the directory for storing sequencing data')
    args = parser.parse_args()
    main(args)