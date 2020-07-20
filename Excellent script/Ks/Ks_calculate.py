# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

import sys
import os
import argparse
fasta_dict = {}
def fasta2dict(fasta_file):
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
               line_list = line.split('>')
               key = line_list[1].split(' ')[0]  #TODO 这里要注意csd文件的header格式！
               fasta_dict[key] = ''
            else:
                fasta_dict[key] += line.replace('\n', '')

    return fasta_dict


# fasta2dict(r'D:\Result\Transit\head5.Sspon.gene.fasta')
# fasta2dict(sys.argv[1])

def get_gene_ortholog_seq(anchor_file):
    if not os.path.exists('01-orthogene_ka_ks_seq'):
        os.mkdir('01-orthogene_ka_ks_seq')
    output_path = '01-orthogene_ka_ks_seq/'
    output = open(output_path + 'Ortholog_gene_kaks_seq.fa', 'w')
    with open(anchor_file) as f:
        for line in f:
            line = line.strip()
            if '#' not in line:
                line_list = line.split('\t')
                for i in line_list[:2]:
                    seq_fa = fasta_dict[i]
                    # print(seq_fa)
                    output.write('>'+i+'\n')
                    output.write(seq_fa+'\n')
    output.close()

# get_gene_ortholog_seq(sys.argv[2])

def run_kaks_calculate():
    if not os.path.exists('02-kaks_result'):
        os.mkdir('02-kaks_result')
    output = open('run_kaks_caculate.sh','w')
    output.write('python /public1/home/stu_panhaoran/scripts/ks_calculate/synonymous_calc.py 01-orthogene_ka_ks_seq/Ortholog_gene_kaks_seq.fa > 02-kaks_result/ka_ks_result.txt ')
    output.close()

# run_kaks_caculate()

def result_merge(result_file):
    path = '02-kaks_result'
    output = open(result_file,'w')
    output.write('name\tdS-yn\tdN-yn\tdS-ng\tdN-ng\n')
    file_name = path +'/'+ 'ka_ks_result.txt'
    with open(file_name) as f:
        for line in f :
            if not line.startswith('name'):
                line = line.replace(',','\t')
                output.write(line)
    output.close()

def main(args):
    if args.fasta and args.anchor :
        fasta_file = args.fasta
        anchor_file = args.anchor
        fasta2dict(fasta_file)
        get_gene_ortholog_seq(anchor_file)
        run_kaks_calculate()
    if args.output:
        result_file = args.output
        result_merge(result_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='Ks_calculate.py',
        description='''Calculate the synonymous mutation rate/non-synonymous mutation rate of two species based on jcvi results''')
    parser.add_argument('-f', '--fasta',help='Input the combined cds.fasta file of the two species')
    parser.add_argument('-a', '--anchor',help='Input the anchor file obtained by jcvi')
    parser.add_argument('-o', '--output',help='Please name the final result')
    args = parser.parse_args()
    main(args)