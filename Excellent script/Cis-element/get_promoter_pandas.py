# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

#不建议使用此方法，因为自己整理的gff3文件末尾容易有空格符，导入pandas时会报错，可使用普通字符脚本进行处理！！
#如果是JGI数据库下载的gff文件可以使用！
import pandas as pd
import argparse


#获得目的基因的基因ID
gene_list = []
def gain_genome_wide_geneID(gene_file):
    with open(gene_file) as f:
        for line in f:
            line = line.rstrip()
            line = line[0:]      # TODO 根据不同物种的geneID进行修改
            gene_list.append(line)
    return gene_list
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
def Operate_on_gff_files(gff_file,skip_rows=1):
    gff = pd.read_table(gff_file,sep='\t',header=None,skiprows=skip_rows)
    gff.rename(
            columns={
                0:'ChrID',1:'Source',2:'Type',3:'Start',4:'End',5:'Score',6:'Strand',7:'Phase',8:'Attributes'},inplace=True)
    cds_gff = gff[gff['Type'] == 'CDS']
    for index,row in cds_gff.iterrows():
        Chr = cds_gff.loc[index]['ChrID']
        strand = cds_gff.loc[index]['Strand']
        if strand == '+':
            start_position = str(cds_gff.loc[index]['Start'])
        else:
            start_position = str(cds_gff.loc[index]['End'])
        # print(start_position)
        CDS_ID = cds_gff.loc[index]['Attributes'].split(';')[1].split('=')[1]     # TODO 根据不同物种的gff文件进行修改
        if CDS_ID not in gff_dict:
            gff_dict[CDS_ID] = Chr + '_' + start_position + '_' + strand

    return gff_dict


#最后的处理及结果文件的输出
def Final_processing_and_output_of_result_files(output_file):
    result = open(output_file, 'w')
    for i in gene_list:
        gene_ID = i.strip()
        if gene_ID != '':
            tmp = gff_dict[gene_ID].split('_')
            if tmp[2] == '+':
                seq = genome_dict[tmp[0]][int(tmp[1])-1:int(tmp[1])-2000-1:-1]
                seq_list = list(seq)
                seq_list.reverse()
                seq = ''.join(seq_list)
                result.write('>' + gene_ID + ' +' + '\n' + seq + '\n')
            else:
                result.write('>' + gene_ID + ' -' + '\n' + genome_dict[tmp[0]][int(tmp[1]):int(tmp[1]) + 2000] + '\n')
    result.close()

def main(args):
    gene_file = args.geneID
    genome_fasta = args.genome
    gff_file = args.gff
    skip_rows = args.skip
    output_file = args.out
    gain_genome_wide_geneID(gene_file)
    Save_genome_files_as_a_dict(genome_fasta)
    Operate_on_gff_files(gff_file,skip_rows)
    Final_processing_and_output_of_result_files(output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='gain_promoter_region.py',
        description='''This script is used to get the genome-wide promoter sequence (up 2000bp) of each gene''')
    parser.add_argument('-i','--geneID',required=True,help='gene_file')
    parser.add_argument('-g','--genome',required=True,help='genome.fasta')
    parser.add_argument('-f','--gff',required=True,help='gff_file')
    parser.add_argument('-s','--skip',type = int,help='gff_file_skiprows')
    parser.add_argument('-o','--out',required=True,help='output_file')
    args = parser.parse_args()
    main(args)
