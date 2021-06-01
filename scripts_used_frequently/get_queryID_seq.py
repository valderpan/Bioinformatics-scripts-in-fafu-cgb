#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/5/31

'''
%prog <fasta> <geneID file> <query geneID fasta>

Find the fasta file of the target gene based on the target gene ID
'''

from Bio import SeqIO
import sys
import SetLog

def parse_format(fasta_file):
    Dict = {}
    for seq in SeqIO.parse(fasta_file,'fasta'):
        Dict[seq.id] = seq.seq
    return Dict

def read_geneID(geneID_file):
    geneid = []
    with open(geneID_file,'r') as f:
        lines = (line.strip() for line in f)
        for line in lines:
            geneid.append(line)
    query_IDnumber = len(geneid)
    return geneid,query_IDnumber

def match_query_fasta(fasta_Dict,geneID_list,IDnumber,output):
    matched_num = 0
    for key in fasta_Dict.keys():
        if key in geneID_list:
            matched_num +=1
    if matched_num != IDnumber:
        SetLog.warning_out('Matcher number != query number !!')
        with open(output,'w') as f1:
            for key in fasta_Dict.keys():
                if key in geneID_list:
                    f1.write('>'+key+'\n')
                    f1.write(str(fasta_Dict[key])+'\n')
    else:
        SetLog.info_out('Matcher number == query number |[check ok]')
        with open(output,'w') as f1:
            for key in fasta_Dict.keys():
                if key in geneID_list:
                    f1.write('>'+key+'\n')
                    f1.write(str(fasta_Dict[key])+'\n')


if __name__ == "__main__":
    from optparse import OptionParser
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args)==3:
        fasta_file = args[0]
        geneID_file = args[1]
        output_file = args[2]
        fasta_Dict = parse_format(fasta_file)
        geneID_list,query_IDnumber = read_geneID(geneID_file)
        match_query_fasta(fasta_Dict,geneID_list,query_IDnumber,output_file)
    else:
        sys.exit(p.print_help())