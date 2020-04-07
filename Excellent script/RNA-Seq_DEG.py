# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

'''
RNA-Seq定量分析
预先创建两个文件夹，一个是01-ref 用来存放参考基因组 ； 另一个是02-reads 用来存放二代测序reads
此脚本默认以用conda安装了hisat2  samtools  featureCounts
'''
import os
import subprocess
import datetime
import re
import argparse

############### MODIFY THE FOLLWINGS PATHS FOR ALL DEPENDENT PROGRAMS ###############
hisat2 = ''
samtools = ''
featurecounts = ''
#####################################################################################


reads_list = []
reads_ID_list = []
def Get_sequencing_data(path):    #path默认是02-reads
    os.chdir(path)
    all_reads = os.listdir('./')
    for read in all_reads:
        read = read.strip()
        read_ID = re.findall('([a-z0-9A-Z\-]+)\_',read)
        reads_list.append(read_ID)
    for ID in reads_list:
        if ID not in reads_ID_list:
            reads_ID_list.append(ID)
    os.chdir('../')
    return reads_ID_list

def run_command(cmd):
    print(cmd)
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
        print("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd,datetime.datetime.now()))

def build_index(genome):
    os.mkdir('03-index')
    cmd1 = 'hisat2-build ' + './01-ref/'+ genome + ' ./03-index/genome'
    run_command(cmd1)

def hisat2_mapping():
    os.mkdir('04-hisat2_results')
    reads_ID_str = [''.join(map(str,sub_list)) for sub_list in reads_ID_list]
    for read in reads_ID_str:
        cmd2 = 'hisat2 ' + '-p 12 ' + '-x ./03-index/genome ' + '-1 '+'./02-reads/'+read+'_R1.fq.gz ' + '-2 '+'./02-reads/' +read+'_R2.fq.gz '+ '-S '+read+'.sam'
        run_command(cmd2)
        cmd3 = 'samtools '+'sort '+'-@ 12 '+ '-o '+read +'.bam '+ read+ '.sam'
        run_command(cmd3)
        cmd4 = 'rm '+read+'.sam'
        run_command(cmd4)
    cmd5 = 'mv ./*.bam ./04-hisat2_results/'
    run_command(cmd5)

def featurecounts_quantitative(gff3_file):
    os.mkdir('05-featurecounts_results')
    cmd6 = 'featureCounts ' + '-T 15 ' + '-p -t gene -g ID ' + '-a '+ './01-ref/'+ gff3_file + ' -o ./05-featurecounts_results/counts_id.txt ' + './04-hisat2_results/*.bam'
    run_command(cmd6)

def main(args):
    reads_path = args.reads
    ref_genome = args.genome
    ref_gff = args.gff

    Get_sequencing_data(reads_path)
    print('Step1 build ref-genome index')
    build_index(ref_genome)

    print('Step2 mapping (This step will take some time)')
    hisat2_mapping()

    print('Step3 reads quantitative')
    featurecounts_quantitative(ref_gff)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='RNA-Seq_DEG.py',
        formatter_class=argparse.RawTextHelpFormatter,
        description='''This script is used to execute RNA-Seq pipeline to find differentially expressed genes''')
    parser.add_argument(
        '-r','--reads',
        required=True,
        help='Specify the directory where the sequencing data is located')
    parser.add_argument(
        '-g','--genome',
        required=True,
        help='Reference genome sequence file')
    parser.add_argument(
        '-f','--gff',
        required=True,
        help='Annotated files for reference genomes')
    args = parser.parse_args()
    main(args)
