# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

import os
from collections import OrderedDict
import argparse
from datetime import datetime
import time
import subprocess


#Modify the path of all dependent programs below
MAFFT = '~/miniconda3/bin/mafft'
RAxML = '~/miniconda3/bin/raxmlHPC-PTHREADS'
TRIMAL = '~/miniconda3/bin/trimal'


SpeciesID = []
SingleCopyGeneID_list = []
SingleCopyGeneRow_list = []
pep_fa_dict = OrderedDict()
seq_aln_trimed_file_dict = OrderedDict()
merge_allSCGene_dict = OrderedDict()
OrthoDict = OrderedDict()


def Get_speciesID(speciesID):
    '''
    Provide abbreviations of species in Orthogroups.tsv
    '''
    with open(speciesID) as f:
        for line in f :
            line = line.strip()
            SpeciesID.append(line)
    return SpeciesID


def run_command(cmd):
    '''
    Run shell command
    '''
    print(cmd)
    return_code = subprocess.call(cmd,shell=True)
    if return_code != 0 :
        print('ERROR:[{2}] Return code {0} when running the following command : {1} '.format(return_code,cmd,datetime.now()))


def Get_SingleCopyGeneID(SingleOrthoGroup,OrthoGroups):
    '''
    Process Orthogroups.tsv according to Orthogroups_SingleCopyOrthologues.txt
    '''
    with open(SingleOrthoGroup) as f :
        for line in f:
            line = line.strip()
            SingleCopyGeneID_list.append(line)
    with open(OrthoGroups) as f :
        for line in f :
            line = line.strip()
            if line.startswith('Ortho'):
                continue
            else:
                line = line.split()
                OrthoDict[line[0]] = '='.join(line[1:])
    for i in SingleCopyGeneID_list:
        if i in OrthoDict:
            SingleCopyGeneRow_list.append(i + ':' + OrthoDict[i])
    return SingleCopyGeneRow_list


def Get_SingleCopyGeneSeq(SingleCopyGeneRow,allseq):
    '''
    Get the sequence of a single copy orthologous gene
    '''
    with open(allseq) as f:
        for line in f:
            if line.startswith(">"):
                seq_header = line.strip('>').split()[0]
                pep_fa_dict[seq_header] = ''
            else:
                pep_fa_dict[seq_header] += line
    os.mkdir('Singlecopygene')
    os.chdir('Singlecopygene')
    for line in SingleCopyGeneRow:
         line = line.strip().split(':')
         seq_file = line[0] + '.fasta'
         gene_ids = line[1].split('=')
         with open(seq_file,'w') as f :
             for i in gene_ids:
                 if i in pep_fa_dict:
                     for j in SpeciesID :
                         if i.startswith(j):
                             f.write('>'+j+'\n'+'\n'.join(pep_fa_dict[i])+'\n')
         f.close()
    os.chdir('../')



def MSA_SingleCopyGene(SingleCopyGeneRow):
    '''
    Run MAFFT and Trimal
    '''
    os.mkdir("SingleGene_MSA")
    os.chdir("SingleGene_MSA")
    for line in SingleCopyGeneRow:
        line = line.strip().split(':')
        seq_file = line[0] + '.fasta'
        seq_aln_file = line[0] + '.aln.fasta'
        seq_aln_trimed_file = line[0] + "_aln_trimed.fas"
        cmd1 = MAFFT + ' ../Singlecopygene/' + seq_file + ' >' + seq_aln_file
        run_command(cmd1)
        cmd2 = TRIMAL + ' -in ' + seq_aln_file + ' -out ' + seq_aln_trimed_file + ' -automated1'
        run_command(cmd2)


def merge_SingleCopyGene(SingleCopyGeneRow):
    '''
    Merge every aligned single-copy gene into one super-gene matirx
    '''
    for line in SingleCopyGeneRow:
        line = line.strip().split(':')
        seq_aln_file = line[0] + '.aln.fasta'
        seq_aln_trimed_file = line[0] + "_aln_trimed.fas"
        with open(seq_aln_trimed_file) as f:
            for i in f:
                if i.startswith('>'):
                    seqid = i.strip(">").split()[0] + "." + line[0]
                    seq_aln_trimed_file_dict[seqid] = ''
                else:
                    seq_aln_trimed_file_dict[seqid] += i
            for ids in SpeciesID:
                merge_allSCGene_dict[ids] = []
                for key, value in seq_aln_trimed_file_dict.items():
                    if key.startswith(ids):
                        merge_allSCGene_dict[ids].append("".join(value))

    with open("merged_allSingleGenes.fas", "w") as fh:
        for key, value in merge_allSCGene_dict.items():
            fh.write(">" + key + "\n" + "".join(value) + "\n")


def MLtree_Concatenation():
    '''
    Construct the ML species tree with the concatenate methold
    '''
    os.mkdir("Concatenation")
    os.chdir("Concatenation")
    cmd = RAxML + ' -T ' + thread + ' -f a -N ' + nb + ' -m ' + model + ' -x 123456 -p 123456 -s ../SingleGene_MSA/merged_allSingleGenes.fas -n concatenation_out.nwk'
    run_command(cmd)
    os.chdir("../")


def main(args):
    reader1 = args.input1
    reader2 = args.input2
    reader3 = args.input3
    reader4 = args.input4
    global thread, nb, model
    thread = args.thread
    nb = args.bootstrap
    model = args.model

    Get_speciesID(reader1)
    print("Step1: Get each single-copy gene id and its sequences.\n")
    single_gene_ids = Get_SingleCopyGeneID(reader2, reader3)
    Get_SingleCopyGeneSeq(single_gene_ids, reader4)

    print("Step2: Conduct multiple sequences alignment for each single-copy gene.\n")
    MSA_SingleCopyGene(single_gene_ids)

    print("Step3: Merge all aligned single-copy gene into one file.\n")
    merge_SingleCopyGene(single_gene_ids)
    os.chdir("../")

    print("Step4: Construct the ML species tree with the concatenate method.\n")
    MLtree_Concatenation()

    print("Congratulations! All tasks were finished!\n")

if __name__ == '__main__':
    start = time.perf_counter()
    parser = argparse.ArgumentParser(
        prog= 'Orthofinder2SpeciesTree',
        formatter_class= argparse.RawTextHelpFormatter,
        description='''
        %(prog)s <SpeciesID prefix> <SingleCopyOrtho> <Orthogroups> <protein file> [thread] [bootstrap] [model]

        Use this script to get phylogenetic tree files from OrthoFinder results files
        ''')
    parser.add_argument(
        '-in1', '--input1',
        required=True,
        help="offer the prefix of all abbreviated species id ")
    parser.add_argument(
        '-in2', '--input2',
        required=True,
        help="offer the Single-copy Orthogroups file, SingleCopyOrthogroups.txt")
    parser.add_argument(
        '-in3', '--input3',
        required=True,
        help="offer the all Orthogroups file, Orthogroups.tsv")
    parser.add_argument(
        '-in4', '--input4',
        required=True,
        help="offer all species protein sequences")
    parser.add_argument(
        '-t', '--thread',
        default="10",
        help="set the number of thread, default=10")
    parser.add_argument(
        '-nb', '--bootstrap',
        default="100",
        help="set the number of bootstrap, default=100")
    parser.add_argument(
        '-m', '--model',
        default="PROTGAMMAJTT",
        help="set the model of amino acid substitution, default=PROTGAMMAJTT")
    args = parser.parse_args()
    main(args)
    end = time.perf_counter()
    print("All tasks used time: %ss" % (end - start))


