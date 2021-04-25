#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/4/24

'''
make species' cds phylip base SingCopyGene
'''

from path import Path
import sys
import pandas as pd
import re
from Bio import SeqIO
import os
import argparse
sys.path.append('D:/pycharm/panhr/FAFU_CGB/Genomic_exercise')
import SetLog as Log


__author__ = 'Haoran Pan'
__mail__ = 'haoran_pan@qq.com'
__date__ = '20210425'
__version__ = 'v1.0'


def get_SCG(pep_msa_path):
    OG2geneID = {}
    msa_files = [file for file in Path(pep_msa_path).files() if file.endswith('msa')]
    for msa in msa_files:
        OG = re.findall('(OG[0-9]+)\.msa',msa.basename())[0]
        geneID = []
        with open(msa,'r') as f:
            lines = (line.strip() for line in f)
            for line in lines:
                if line.startswith('>'):
                    geneID.append(line[1:])
        OG2geneID[OG] = geneID
        # break

    return OG2geneID,msa_files


def cds2dict(cds_file_list):
    Allcds_ID2seq = {}
    cdsID = []
    for cds in cds_file_list:
        for seq in SeqIO.parse(cds, 'fasta'):
            Allcds_ID2seq[seq.id] = seq.seq
            cdsID.append(seq.id)
    return Allcds_ID2seq,cdsID


def extract_CDSseq(OG2geneID,cds_dict,cdsID):
    '''
    This function finds the corresponding ID in the cds sequence according to pepID
    '''
    if len([dir for dir in Path('./').dirs() if dir.basename()=='cds_SCG']) == 0:
        os.mkdir('../cds_SCG')

    for key in OG2geneID.keys():
        with open('./cds_SCG/'+key+'.fa','w') as f:
            for ID in OG2geneID[key]:  #ID:pep_ID   cdsID:cds_ID
                if ID in cdsID:
                    f.write('>' + ID + '\n')
                    f.write(str(cds_dict[ID]) + '\n')
                elif ID.startswith('LOC_Os'):           # Osativa
                    if ''.join(ID.split('.')[0:2]) in cdsID:
                        f.write('>' + ID + '\n')
                        f.write(str(cds_dict[''.join(ID.split('.'))]) + '\n')
                elif ID.startswith('Chr'):       # Osativa
                    if re.sub('mRNA','gene',ID) in cdsID:
                        f.write('>' + ID + '\n')
                        f.write(str(cds_dict[re.sub('mRNA','gene',ID)]) + '\n')
                elif ID.startswith('Sobic'):        #Sbicolor
                    if '.'.join(ID.split('.')[0:3]) in cdsID:
                        f.write('>' + ID + '\n')
                        f.write(str(cds_dict['.'.join(ID.split('.')[0:3])]) + '\n')
                elif ID.startswith('Zm'):       #Zea Mays
                    if ID.split('_')[0] in cdsID:
                        f.write('>' + ID + '\n')
                        f.write(str(cds_dict[ID.split('_')[0]]) + '\n')


def debug_cdsfa(msafiles,species_num):
    cdsfafiles = [file for file in Path('./cds_SCG/').files() if file.endswith('.fa')]
    if not len(msafiles) == len(cdsfafiles):
        Log.errot_out('The number of generated cds fasta files is inconsistent with the number of single-copy protein files')
    else:
        Log.debug_out('Correct ! The number of generated cds fasta files is the same as the number of single-copy protein files !')
    for cds in cdsfafiles:
        num = 0
        with open(cds,'r') as f:
            lines = (line.strip() for line in f)
            for line in lines:
                if line.startswith('>'):
                    num +=1
            if not num == species_num:
                Log.error_out('The number of sequences in the {} file is not {} !'.format(cds,species_num))
            else:
                continue


def main(args):
    pepmsa_path = args.pepmsa
    cdsseq = args.cdsfiles
    species_num = args.number
    OG2ID,msafiles = get_SCG(pepmsa_path)
    cds_dict,cdsID = cds2dict(cdsseq)
    extract_CDSseq(OG2ID,cds_dict,cdsID)
    debug_cdsfa(msafiles,species_num)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''make species' cds phylip base SingCopyGene''',
        usage="python {} -a -o1 ".format(sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__, __mail__, __date__,
                                                                            __version__))
    parser.add_argument('--pepmsa', required=True, help='Input pep msa path')
    parser.add_argument('-c', '--cdsfiles',required=True,nargs='+', help='Input species whole genome cds fasta')
    parser.add_argument('-n', '--number', type=int, help='Input the number of specie')

    args = parser.parse_args()
    main(args)