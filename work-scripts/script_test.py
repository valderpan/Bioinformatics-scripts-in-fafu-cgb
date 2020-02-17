# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# 下面是师兄的fasta_to_dictionary.py
# def fasta_to_dict():
#     import sys
#     fasta_seq = open(sys.argv[2])
#     fasta_line = fasta_seq.readline()
#     fasta_dict = {}
#     while fasta_line:
#         if fasta_line.startswith('>'):
#             value_seq = fasta_line.rstrip()
#             fasta_line_list = fasta_line.split(">")
#             key = fasta_line_list[1].rstrip()
#             fasta_line = fasta_seq.readline()
#             while fasta_seq :
#                 if not fasta_line.startswith('>'):
#                     value_seq = value_seq + '\n' + fasta_line.rstrip()
#                     fasta_line = fasta_seq.readline()
#                 else :
#                     fasta_dict[key] = value_seq
#                     break
#         else :
#             fasta_line = fasta_seq.readline()
#     return fasta_dict
#
# if __name__ == '__main__':
#     fasta_to_dict()

#下面是按我自己的习惯做了一遍
# import os
# os.chdir(r'D:\Result\test')
# fa_dict = {}
# # def fasta_to_dict():
# with open('reads_2.fa') as f :
#     fasta_line = f.readline()
#     # fa_dict = {}
#     while fasta_line: #第一行进入fasta_line
#         if fasta_line.startswith('>'):
#             value_seq = fasta_line.rstrip()
#             key = fasta_line[1:].rstrip()
#             fasta_line = f.readline() #第二行进入fasta_line
#             while fasta_line:
#                 if not fasta_line.startswith('>'):
#                     value_seq = value_seq+ '\n'+fasta_line.rstrip()
#                     fasta_line = f.readline()
#                 else:
#                     fa_dict[key] = value_seq
#                     # print(fa_dict[key])
#                     break
#         else:
#             fasta_line=f.readline()
#             # print(fasta_line)
#             # break
#     # print(fa_dict['EVM0010950.1'])
#     print(fa_dict.keys())

# # fasta_to_dict()

# 下面是按照我自己的习惯测试的师兄的01.creat_ortholog.py
import os
# import fasta_to_dictionary
os.chdir(r'D:\Result')
with open('yellowhorn.Xs.anchors') as f:
    dict_fa = {}
    for line  in f :
        line1 = line.strip()
        if not '#' in line1:
            line1_list = line1.split('\t')
            for line2 in line1_list[:2]:
                seq_fa = dict_fa[line2.rstrip()]
                print(seq_fa)