# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

import argparse

# Description:This script is used to get the genome-wide promoter sequence of each gene


#获得目的基因的基因ID
gene_list = []
def gain_genome_wide_geneID(gene_file):
    with open(gene_file) as f:
        for line in f:
            line = line.rstrip()
            line = line[0:]   #TODO 根据物种特定的基因ID进行更改
            gene_list.append(line)

#将基因组文件储存为字典
genome_dict = {}
def Save_genome_files_as_a_dict(genome_fasta):
    with open(genome_fasta) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                key = line[1:]
                genome_dict[key] = line
            else:
                line = line.upper()
                genome_dict[key] += line.replace('\n','')
    return genome_dict

#对gff注释文件进行操作
gff_dict = {}
def Operate_on_gff_files(gff_file):
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
                    CDS_id = '-'.join(line_list[8].split(';')[0].split('=')[-1].split('-')[0:2])  #TODO 根据物种gff文件中的CDSID进行更改
                    # print(CDS_id)
                    if CDS_id not in gff_dict:
                        gff_dict[CDS_id] = Chr + '_' + start_position + '_'+strand
    return gff_dict


#最后的处理及结果文件的输出
def Final_processing_and_output_of_result_files(output_file):
    result = open(output_file, 'w')
    for i in gene_list:
        gene_ID = i.rstrip()
        if gene_ID != '':
            tmp = gff_dict[gene_ID].split('_')
            if tmp[2] == '+':
                seq = genome_dict[tmp[0]][int(tmp[1])-1:int(tmp[1])-2000-1:-1]
                seq_list = list(seq)
                seq_list.reverse()
                seq = ''.join(seq_list)
                result.write('>' + gene_ID + ' +' + '\n'+ seq + '\n')
            else:
                seq = genome_dict[tmp[0]][int(tmp[1]):int(tmp[1])+2000]
                seq_list = list(seq)
                seq_list.reverse()
                seq = ''.join(seq_list)
                result.write('>' + gene_ID + ' -' + '\n'+ seq + '\n')

def main(args):
    gene_file = args.geneID
    genome_fasta = args.genome
    gff_file = args.gff
    output_file = args.out
    gain_genome_wide_geneID(gene_file)
    Save_genome_files_as_a_dict(genome_fasta)
    Operate_on_gff_files(gff_file)
    Final_processing_and_output_of_result_files(output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='gain_promoter_region.py',
        description='''This script is used to get the genome-wide promoter sequence of each gene''')
    parser.add_argument('-i','--geneID',required=True,help='gene_file')
    parser.add_argument('-g','--genome',required=True,help='genome.fasta')
    parser.add_argument('-f','--gff',required=True,help='gff_file')
    parser.add_argument('-o','--out',required=True,help='output_file')
    args = parser.parse_args()
    main(args)
