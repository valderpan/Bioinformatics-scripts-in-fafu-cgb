# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

#此脚本用于根据orthofinder中的单拷贝同源基因ID将单拷贝同源蛋白序列提取出来，用于多序列比对。

import os
import argparse

SingleCopyOrtho_ID = []
SingleCopyOrtho_group = []
fasta_dict = {}



def Get_SCO_group(Orthogroups,SingleCopyOrthologues):
    with open(SingleCopyOrthologues)as f:
        for ID in f:
            ID = ID.strip()
            SingleCopyOrtho_ID.append(ID)
    with open(Orthogroups) as f:
        for line in f:
            line = line.strip()
            for id in SingleCopyOrtho_ID:
                if id in line:
                    SingleCopyOrtho_group.append(line)
    return SingleCopyOrtho_group


def Get_Single_copy_ortholog_prot(all_pep):
    os.mkdir("SingleCopyOrthogroupseq")
    os.chdir("SingleCopyOrthogroupseq")
    with open('../'+all_pep) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                key = line[1:]
                fasta_dict[key] = ''
            else:
                fasta_dict[key] += line.replace('\n','')

        for line in SingleCopyOrtho_group:
            line =line.strip()
            line = line.split('\t')
            seq_file = open(line[0]+'.fasta','w')
            for j in line[1:]:
                seq_file.write('>'+j + '\n'+fasta_dict[j] + '\n')


def main(args):
    Orthogroups = args.groups
    SingleCopyOrthologues = args.ID
    all_pep = args.all
    Get_SCO_group(Orthogroups,SingleCopyOrthologues)
    Get_Single_copy_ortholog_prot(all_pep)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='gain_SCO_prot2MAFFT.py',
        formatter_class=argparse.RawTextHelpFormatter,
        description='''This script is used to extract single copy orthologous protein sequences from orthofinfer for next multiple sequence alignment.''')
    parser.add_argument('-g','--groups',required=True,help="Enter orthofinder results : Orthogroups.tsv")
    parser.add_argument('-i','--ID',required=True,help="Enter orthofinder results : Orthogroups_SingleCopyOrthologues.txt")
    parser.add_argument('-a','--all',required=True,help="Please combine the protein files of all species into one file and input")
    args = parser.parse_args()
    main(args)