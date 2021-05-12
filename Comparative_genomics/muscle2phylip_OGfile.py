#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/5/12

import muscle2phylip_path
import sys
from Bio import SeqIO
import re
import argparse
sys.path.append('/public1/home/stu_panhaoran/biosoft/pyhr/scripts')
import SetLog as log


__author__ = 'Haoran Pan'
__mail__ = 'haoran_pan@qq.com'
__date__ = '20210512'
__version__ = 'v1.2'


def msa2phylip_OG(msa_files,species_list,output_path):
    for file in msa_files:
        file_prefix = re.findall('(OG[0-9]+)\.msa',file)[0]
        gene2seq = {}
        for seq in SeqIO.parse(file,'fasta'):
            gene2seq[seq.id] = seq.seq
        species_length = [len(gene2seq[key]) for key in gene2seq.keys()]
        length_uniq = list(set(species_length))
        if len(length_uniq) == 1:
            with open('{}/{}.phylip'.format(output_path, file_prefix), 'w') as fw:
                fw.write('{}\t{}\n'.format(len(species_list), length_uniq[0]))
                for key in gene2seq.keys():
                    for specie in species_list:
                        if key.startswith(specie):
                            fw.write('{}    {}\n'.format(specie,gene2seq[key]))

        else:
            log.error_out('Species msa length is inconsistent')

def main(args):
    msapath = args.path
    specieslist = args.species
    outputpath = args.outputpath
    L = muscle2phylip_path.get_msa(msapath)
    msa2phylip_OG(L,specieslist,outputpath)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Convert the result of each multiple sequence alignment to phylip format''',
        usage="python {} -p msa_path -o msa_phylip_path -s ID1 ID2 ID3..\n\nNote:ID is the ID of the gene in the genome of each species!!!\nGenerated file: SingleCopy.pep.phylip".format(
            sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__, __mail__, __date__,
                                                                            __version__))
    parser.add_argument('-p', '--path', required=True, help='Input the path where the msa file is located')
    parser.add_argument('-o', '--outputpath', required=True,
                        help='Specify output file directory')
    parser.add_argument('-s', '--species', required=True, nargs='+',
                        help='Input the gene ID prefix for each species genome')
    args = parser.parse_args()
    main(args)
