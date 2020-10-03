#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# @Time : 2020/10/1 9:54

'''
Usage:python families_count2venn_v3.py -O C:\\Users\\dell\\Desktop\\Nepal_project\\20200525_venn\\0930\\Orthogroups.tsv -k1 Misin -o1 Misin.txt -k2 XSsp -o2 Nepal.txt -k3 LOC_Os -o3 Os.txt -k4 Sspon -o4 Sspon.txt -k5 Sobic -o5 Sb.txt
'''
import argparse
import re

__author__ = 'Haoran Pan'
__mail__ = 'haoran_pan@qq.com'
__version__ = 'v3.0'

def Orthogroups_count2venn_2(Orthogroups_tsv,key1,key2,output1,output2):
    with open(Orthogroups_tsv) as f:
        key1_list = [];key2_list = []
        for line in f:
            line = line.strip()
            if line.startswith('Orthogroup'):
                continue
            else:
                line_list = line.split()
                for i in line_list:
                    if key1 in i:
                        if line_list[0] not in key1_list:
                            key1_list.append(line_list[0])
                for i in line_list:
                    if key2 in i:
                        if line_list[0] not in key2_list:
                            key2_list.append(line_list[0])
        key1_output = open(output1,'w')
        key2_output = open(output2,'w')
        for i in key1_list:
            key1_output.write(i+'\n')
        for i in key2_list:
            key2_output.write(i + '\n')
        key1_output.close()
        key2_output.close()

def Orthogroups_count2venn_3(Orthogroups_tsv,key1,key2,key3,output1,output2,output3):
    with open(Orthogroups_tsv) as f:
        key1_list = [];key2_list = [];key3_list = []
        for line in f:
            line = line.strip()
            if line.startswith('Orthogroup'):
                continue
            else:
                line_list = line.split()
                for i in line_list:
                    if key1 in i:
                        if line_list[0] not in key1_list:
                            key1_list.append(line_list[0])
                for i in line_list:
                    if key2 in i:
                        if line_list[0] not in key2_list:
                            key2_list.append(line_list[0])
                for i in line_list:
                    if key3 in i:
                        if line_list[0] not in key3_list:
                            key3_list.append(line_list[0])
        key1_output = open(output1,'w')
        key2_output = open(output2,'w')
        key3_output = open(output3,'w')
        for i in key1_list:
            key1_output.write(i+'\n')
        for i in key2_list:
            key2_output.write(i + '\n')
        for i in key3_list:
            key3_output.write(i + '\n')
        key1_output.close()
        key2_output.close()
        key3_output.close()


def Orthogroups_count2venn_4(Orthogroups_tsv,key1,key2,key3,key4,output1,output2,output3,output4):
    with open(Orthogroups_tsv) as f:
        key1_list = [];key2_list = [];key3_list = [];key4_list = []
        for line in f:
            line = line.strip()
            if line.startswith('Orthogroup'):
                continue
            else:
                line_list = line.split()
                for i in line_list:
                    if key1 in i:
                        if line_list[0] not in key1_list:
                            key1_list.append(line_list[0])
                for i in line_list:
                    if key2 in i:
                        if line_list[0] not in key2_list:
                            key2_list.append(line_list[0])
                for i in line_list:
                    if key3 in i:
                        if line_list[0] not in key3_list:
                            key3_list.append(line_list[0])
                for i in line_list:
                    if key4 in i:
                        if line_list[0] not in key4_list:
                            key4_list.append(line_list[0])
        key1_output = open(output1,'w')
        key2_output = open(output2,'w')
        key3_output = open(output3,'w')
        key4_output = open(output4,'w')
        for i in key1_list:
            key1_output.write(i+'\n')
        for i in key2_list:
            key2_output.write(i + '\n')
        for i in key3_list:
            key3_output.write(i + '\n')
        for i in key4_list:
            key4_output.write(i + '\n')
        key1_output.close()
        key2_output.close()
        key3_output.close()
        key4_output.close()


def Orthogroups_count2venn_5(Orthogroups_tsv,key1,key2,key3,key4,key5,output1,output2,output3,output4,output5):
    with open(Orthogroups_tsv) as f:
        key1_list = [];key2_list = [];key3_list = [];key4_list = [];key5_list = []
        for line in f:
            line = line.strip()
            if line.startswith('Orthogroup'):
                continue
            else:
                line_list = line.split()
                for i in line_list:
                    if key1 in i:
                        if line_list[0] not in key1_list:
                            key1_list.append(line_list[0])
                for i in line_list:
                    if key2 in i:
                        if line_list[0] not in key2_list:
                            key2_list.append(line_list[0])
                for i in line_list:
                    if key3 in i:
                        if line_list[0] not in key3_list:
                            key3_list.append(line_list[0])
                for i in line_list:
                    if key4 in i:
                        if line_list[0] not in key4_list:
                            key4_list.append(line_list[0])
                for i in line_list:
                    if key5 in i:
                        if line_list[0] not in key5_list:
                            key5_list.append(line_list[0])

        key1_output = open(output1,'w')
        key2_output = open(output2,'w')
        key3_output = open(output3,'w')
        key4_output = open(output4,'w')
        key5_output = open(output5,'w')
        for i in key1_list:
            key1_output.write(i+'\n')
        for i in key2_list:
            key2_output.write(i + '\n')
        for i in key3_list:
            key3_output.write(i + '\n')
        for i in key4_list:
            key4_output.write(i + '\n')
        for i in key5_list:
            key5_output.write(i + '\n')
        key1_output.close()
        key2_output.close()
        key3_output.close()
        key4_output.close()
        key5_output.close()

def Orthogroups_count2venn_6(Orthogroups_tsv,key1,key2,key3,key4,key5,key6,output1,output2,output3,output4,output5,output6):
    with open(Orthogroups_tsv) as f:
        key1_list = [];key2_list = [];key3_list = [];key4_list = [];key5_list = [];key6_list = []
        for line in f:
            line = line.strip()
            if line.startswith('Orthogroup'):
                continue
            else:
                line_list = line.split()
                for i in line_list:
                    if key1 in i:
                        if line_list[0] not in key1_list:
                            key1_list.append(line_list[0])
                for i in line_list:
                    if key2 in i:
                        if line_list[0] not in key2_list:
                            key2_list.append(line_list[0])
                for i in line_list:
                    if key3 in i:
                        if line_list[0] not in key3_list:
                            key3_list.append(line_list[0])
                for i in line_list:
                    if key4 in i:
                        if line_list[0] not in key4_list:
                            key4_list.append(line_list[0])
                for i in line_list:
                    if key5 in i:
                        if line_list[0] not in key5_list:
                            key5_list.append(line_list[0])
                for i in line_list:
                    if key6 in i:
                        if line_list[0] not in key6_list:
                            key6_list.append(line_list[0])
        key1_output = open(output1,'w')
        key2_output = open(output2,'w')
        key3_output = open(output3,'w')
        key4_output = open(output4,'w')
        key5_output = open(output5,'w')
        key6_output = open(output6,'w')
        for i in key1_list:
            key1_output.write(i+'\n')
        for i in key2_list:
            key2_output.write(i + '\n')
        for i in key3_list:
            key3_output.write(i + '\n')
        for i in key4_list:
            key4_output.write(i + '\n')
        for i in key5_list:
            key5_output.write(i + '\n')
        for i in key6_list:
            key6_output.write(i + '\n')
        key1_output.close()
        key2_output.close()
        key3_output.close()
        key4_output.close()
        key5_output.close()
        key6_output.close()

def main(args):
    Species = args.species
    Orthogroups_tsv = args.Orthogroups
    key1 = args.key1
    key2 = args.key2
    output1 = args.output1
    output2 = args.output2
    if Species == 2:
        Orthogroups_count2venn_2(Orthogroups_tsv,key1,key2,output1,output2)
    elif Species == 3:
        key3 = args.key3
        output3 = args.output3
        Orthogroups_count2venn_3(Orthogroups_tsv,key1,key2,key3,output1,output2,output3)
    elif Species == 4:
        key3 = args.key3
        key4 = args.key4
        output3 = args.output3
        output4 = args.output4
        Orthogroups_count2venn_4(Orthogroups_tsv,key1,key2,key3,key4,output1,output2,output3,output4)
    elif Species == 5:
        key3 = args.key3
        key4 = args.key4
        key5 = args.key5
        output3 = args.output3
        output4 = args.output4
        output5 = args.output5
        Orthogroups_count2venn_5(Orthogroups_tsv,key1,key2,key3,key4,key5,output1,output2,output3,output4,output5)
    elif Species == 6:
        key3 = args.key3
        key4 = args.key4
        key5 = args.key5
        key6 = args.key6
        output3 = args.output3
        output4 = args.output4
        output5 = args.output5
        output6 = args.output6
        Orthogroups_count2venn_6(Orthogroups_tsv,key1,key2,key3,key4,key5,key6,output1,output2,output3,output4,output5,output6)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='test.py',
        formatter_class=argparse.RawTextHelpFormatter,
        description='''test script''',
        epilog= 'author:\t{0}\nmail:\t{1}\nversion:\t{2}'.format(__author__,__mail__,__version__))
    parser.add_argument('-s', '--species', required=True,type=int, help='Inter the number of species used to compare')
    parser.add_argument('-O', '--Orthogroups', required=True, help='Input orthofinder output result Orthogroups.tsv')
    parser.add_argument('-k1', '--key1', required=True, help='Inter the gene ID of the first species (shared id)(eg:NP-X:XSsp)')
    parser.add_argument('-k2', '--key2', required=True, help='Inter the gene ID of the second species (shared id)')
    parser.add_argument('-k3', '--key3', help='Inter the gene ID of the third species (shared id)')
    parser.add_argument('-k4', '--key4',  help='Inter the gene ID of the fourth species (shared id)')
    parser.add_argument('-k5', '--key5',  help='Inter the gene ID of the fifth species (shared id)')
    parser.add_argument('-k6', '--key6',  help='Inter the gene ID of the sixth species (shared id)')
    parser.add_argument('-o1', '--output1', required=True, help='Output the result file of the first species (OG number)(.txt)')
    parser.add_argument('-o2', '--output2', required=True, help='Output the result file of the second species (OG number)(.txt)')
    parser.add_argument('-o3', '--output3', help='Output the result file of the third species (OG number)(.txt)')
    parser.add_argument('-o4', '--output4',  help='Output the result file of the fourth species (OG number)(.txt)')
    parser.add_argument('-o5', '--output5',  help='Output the result file of the fifth species (OG number)(.txt)')
    parser.add_argument('-o6', '--output6',  help='Output the result file of the sixth species (OG number)(.txt)')
    args = parser.parse_args()
    main(args)
