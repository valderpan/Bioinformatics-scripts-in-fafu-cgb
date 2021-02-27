#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2020/11/28

import sys
import gzip
import time
sys.path.append('/share/home/stu_panhaoran/scripts')
import Fontcolor
import argparse
from Bio import SeqIO
from SetLog import debug_out,info_out,warning_out,error_out,critical_out

__author__ = 'Haoran Pan'
__mail__ = 'haoran_pan@qq.com'
__date__ = '20210227'
__version__ = 'v1.3'



#--------------step1 计算每个contig的depth---------------------
def depth_count(depthfile,outputfile):
    '''
    计算每个contig的覆盖深度
    :param depthfile: samtools depth -aa -b contig1.scale *.sorted.bam > depthfile
    :param outputfile:
    :return:
    '''
    reads_num = 0
    coverage_depth = 0
    contig = ''
    result = open(outputfile,'w')
    with open(depthfile) as f:
        for line in f:
            line_list = line.split('\t')
            if reads_num == 0:
                coverage_depth += int(line_list[2].rstrip())
                reads_num += 1
                contig = line_list[0].rstrip()
            elif reads_num != 0 and line_list[0] == contig:
                coverage_depth += int(line_list[2].rstrip())
                reads_num += 1
            else:
                depth = float(coverage_depth) / float(reads_num)
                reads_num = 0
                coverage_depth = 0
                output = contig + '\t' + str(depth)
                # print(output.rstrip())
                result.write(output.rstrip()+'\n')
        last_depth = float(coverage_depth)/float(reads_num)
        # print(contig+'\t'+str(last_depth))
        result.write(contig+'\t'+str(last_depth)+'\n')
    result.close()





#--------------step2 计算标准测序深度-------------------------

def calcultate_sample_length(fq_gz_file):
    '''
    计算二代reads的一个sample双端文件的bp长度，用于计算标准深度
    若一个双端测序文件深度不够，可将两个R1合并为mergeR1输入
    :param fq_gz_file:
    :return:
    '''
    start = time.time()
    reads_num = 0
    line_index = 0
    reads_length = 0
    with gzip.open(fq_gz_file) as fq1:
        for line in fq1:
            line_index +=1
            if line_index %4 == 1:
                reads_num += 1
            elif line_index % 4 ==2:
                reads_length += len(line)
    end = time.time()
    # print('File {} reads num are {}'.format(str(fq_gz_file),reads_num))
    Fontcolor.Tips_output('File {} reads num are {}'.format(str(fq_gz_file),reads_num))
    # print('File {} reads length are {}'.format(str(fq_gz_file),reads_length))
    Fontcolor.Tips_output('File {} reads length are {}'.format(str(fq_gz_file),reads_length))
    # print('Time cost: {:.5f}'.format(end-start))
    Fontcolor.Tips_output('Time cost: {:.5f}'.format(end-start))
    return reads_length



def cal_genome_size(genome_fasta):
    '''
    计算三代contig genome的基因组大小
    :param genome_fasta:
    :return:
    '''
    fasta_dict ={}
    for seq in SeqIO.parse(genome_fasta,'fasta'):
        fasta_dict[seq.id] = seq.seq.upper()
    genome_size = 0
    for key in fasta_dict.keys():
        # print(key,fasta_dict[key],sep='\n')
        genome_size += len(fasta_dict[key])
    Fontcolor.Tips_output('Genome size is {} bp'.format(genome_size))
    return genome_size
        

def run_step2(fqfile,genomefile):
    '''
    用二代reads的总长度/genome的大小得出标准测序深度
    :param fqfile:
    :param genomefile:
    :return:
    '''
    fq_length = calcultate_sample_length(fqfile)
    genomesize = cal_genome_size(genomefile)
    std_depth = round((fq_length *2) / genomesize)
    Fontcolor.Tips_output('Standard Sequencing Depth is {}X'.format(std_depth))
    return std_depth


#--------------step3 计算标准测序深度-------------------------

def Fa2dict(fasta_file):
    '''
    将genome.fasta存入字典
    :param fasta_file:
    :return:
    '''
    fasta_dict ={}
    for seq in SeqIO.parse(fasta_file,'fasta'):
        fasta_dict[seq.id] = seq.seq.upper()
    return fasta_dict

        


def ctg_phasing_base_depth(fasta_dict,depthcountfile,std_depth,Ploidy,result):
    '''
    根据标准深度与每个contig的深度得到需要phasing的congtigID，将其序列复制若干次
    :param fasta_dict:
    :param depthcountfile:
    :param std_depth:
    :param result:
    :return:
    '''
    # ctg_file   = fasta_dict
    depth_file = open(depthcountfile,'r')
    std_depth  = int(std_depth)
    ploidy = int(Ploidy)
    phasing_file=open('Phasing_record.txt','w')
    last_output = open(result,'w')
    for dep_line in depth_file:
            dep_line_list=dep_line.rstrip().split('\t')
            raw_depth=float(dep_line_list[1])/float(std_depth)
            depth=int(round(raw_depth))
            if depth > 1 :
                    phasing_file.write(dep_line.rstrip()+'\t'+str(raw_depth)+'\n')
                    if depth - 1 < ploidy:
                        for i in range(depth-1):
                                # print('>'+dep_line_list[0]+'_'+str(i))
                                last_output.write('>'+dep_line_list[0]+'_'+str(i)+'\n')
                                # print(fasta_dict[dep_line_list[0]])
                                last_output.write(str(fasta_dict[dep_line_list[0]])+'\n')
                    else:
                        for i in range(ploidy - 1):
                            last_output.write('>' + dep_line_list[0] + '_' + str(i) + '\n')
                            last_output.write(str(fasta_dict[dep_line_list[0]]) + '\n')
    phasing_file.close()
    last_output.close()
    depth_file.close()


def main(args):
    depthfile = args.depth
    depthfileout = args.depthout
    info_out('-' * 10 + 'step 1 beginning !' + '-' * 10)
    depth_count(depthfile,depthfileout)
    info_out('-' * 10 + 'step 1 completed !' + '-' * 10)

    if not args.stdth:
        fq_file = args.fastq #二代的reads,R1\R2中的一个即可
        genome_file = args.genome #三代的基因组文件
        info_out('-' * 10 + 'step 2 beginning !' + '-' * 10)
        std_depth = run_step2(fq_file,genome_file)
        info_out('-' * 10 + 'step 2 completed !' + '-' * 10)
    else:
        info_out('-' * 10 + 'step 2 beginning !' + '-' * 10)
        std_depth = args.stdth
        Fontcolor.Tips_output('Input std depth is {}X'.format(std_depth))
        info_out('-' * 10 + 'step 2 completed !' + '-' * 10)

    Ploidy = args.ploidy
    output = args.output
    genome_file = args.genome
    info_out('-' * 10 + 'step 3 beginning !' + '-' * 10)
    fasta_dict = Fa2dict(genome_file)
    ctg_phasing_base_depth(fasta_dict,depthfileout,std_depth,Ploidy,output)
    info_out('-' * 10 + 'step 3 completed !' + '-' * 10)


if __name__ == '__main__':
    # depth_count(sys.argv[1],sys.argv[2])
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Polyploid genome phasing:Use Second-generation sequencing reads mapping Third-generation sequencing genome''',
        usage="python {} -d depth -o depthout -f NGS2_fq.gz -g genome -p ploidy -O output\n\t\tor\n\t\tpython {} -d depth -o depthout -s std_depth -g genome -p ploidy -O output".format(sys.argv[0],sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__, __mail__, __date__,
                                                                            __version__))
    parser.add_argument('-d', '--depth', required=True, help='Input samtools depth outputfile')
    parser.add_argument('-o', '--depthout', required=True, help='Output ctg depth result(.tab)')
    parser.add_argument('-f', '--fastq', required=None, help='Input sample file(R1.fq.gz or R2.fq.gz)')
    parser.add_argument('-s', '--stdth', required=None, help='Input the standard depth of known sequencing data, which is equal to the sequencing data/genome size(bp)')
    parser.add_argument('-g', '--genome', required=True, help='Input Thr Sequencing assmbly genome(.fasta)')
    parser.add_argument('-p', '--ploidy', required=True, help='Input Ploidy of species(eg:8)')
    parser.add_argument('-O', '--output', required=True, help='Output the fasta sequence of the contig that only needs to be phased(.fasta)')

    args = parser.parse_args()
    main(args)
