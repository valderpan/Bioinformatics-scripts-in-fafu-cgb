#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2020/12/19


import sys
import Fontcolor as Font

def Fa2dict(fasta_file,output_fasta=0):
    '''
    此函数将">"号后面的所有内容存为key
    :param fasta_file:
    :param output_fasta:
    :return:
    '''
    fa_dict = {}
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                fa_key = line.strip().split('>')[-1]
                fa_dict[fa_key] = ''
            else:
                fa_dict[fa_key] += line.strip().upper().replace('\n', '')
    if output_fasta:
        output=open(output_fasta,'w')
        for key in fa_dict.keys():
            output.write('>'+key+'\n')
            output.write(fa_dict[key]+'\n')
        output.close()
    return fa_dict

if __name__ == '__main__':
    if len(sys.argv) != 3:
        Font.Scripts_tip('Store the fasta in the dictionary and output the sequence in one line')
        print('Usage:')
        print('\tpython {0} <in.fasta> <out.fasta>'.format(sys.argv[0]))
    else:
        Fa2dict(sys.argv[1],sys.argv[2])
