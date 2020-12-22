# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2020/12/19

import sys
sys.path.append('../')
import argparse
import fa2dict


def gain_genome_wide_geneID(gene_file):
    '''
    获得目的基因的基因ID
    根据gene.fasta或cds.fasta的header或者gff3文件获得全基因组基因ID
    :param gene_file:
    :return:
    '''
    with open(gene_file) as f:
        gene_list = [line.strip()[0:] for line in f]
    return gene_list


def Operate_on_gff_files(gff_file):
    '''
    对gff注释文件进行操作,最终存储为 GeneID:Chr01A(chr)@1000(start)1@+(strand)
    :param gff_file:
    :return:
    '''
    gff_dict = {}
    with open(gff_file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('#'):
                continue
            else:
                line_list = line.split()
                if len(line_list)>2 and line_list[2] == 'CDS':
                    Chr = line_list[0]
                    strand = line_list[6]
                    if strand == '+':
                        start_position = line_list[3]
                    else:
                        start_position = line_list[4]
                    # TODO 根据物种gff文件中的CDSID进行更改  下面输出的CDS_id为Sspon.02G0010910-3D(geneid)
                    CDS_id = '-'.join(line_list[8].split(';')[0].split('=')[-1].split('-')[0:2])
                    # print(CDS_id)
                    if CDS_id not in gff_dict:
                        gff_dict[CDS_id] = Chr + '@' + start_position + '@'+strand
    return gff_dict


def Final_processing_and_output_of_result_files(genelist,genomedict,gffdict,output_file):
    '''
    最后的处理及结果文件的输出
    :param genelist:
    :param genomedict:
    :param gffdict:
    :param output_file:
    :return:
    '''
    result = open(output_file, 'w')
    for i in genelist:
        gene_ID = i.rstrip()
        # print(gene_ID)
        if gene_ID != '':
            tmp = gffdict[gene_ID].split('@')
            if tmp[2] == '+':
                seq = genomedict[tmp[0]][int(tmp[1])-1:int(tmp[1])-2000-1:-1]
                seq_list = list(seq)
                seq_list.reverse()
                seq = ''.join(seq_list)
                result.write('>' + gene_ID + ' +' + '\n'+ seq + '\n')
            else:
                seq = genomedict[tmp[0]][int(tmp[1]):int(tmp[1])+2000]
                seq_list = list(seq)
                seq_list.reverse()
                seq = ''.join(seq_list)
                result.write('>' + gene_ID + ' -' + '\n'+ seq + '\n')
    result.close()


def main(args):
    gene_file = args.geneID
    genome_fasta = args.genome
    gff_file = args.gff
    output_file = args.out
    genelist = gain_genome_wide_geneID(gene_file)
    genomefasta = fa2dict.Fa2dict(genome_fasta)
    gffdict = Operate_on_gff_files(gff_file)
    Final_processing_and_output_of_result_files(genelist,genomefasta,gffdict,output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''This script is used to get the genome-wide promoter sequence(2000bp) of each gene''')
    parser.add_argument('-i','--geneID',required=True,help='gene_file')
    parser.add_argument('-g','--genome',required=True,help='genome.fasta')
    parser.add_argument('-f','--gff',required=True,help='gff_file')
    parser.add_argument('-o','--out',required=True,help='output_file')
    args = parser.parse_args()
    main(args)
