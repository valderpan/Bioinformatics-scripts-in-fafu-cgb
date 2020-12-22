# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-


import os
import sys
sys.path.append('../')
import fa2dict
import Fontcolor as Font


def count_base_number(fa_dict):
    print('GeneID','A_ratio','T_ratio','C_ratio','G_ratio',sep='\t')
    for i in fa_dict.keys():
        A_count = fa_dict[i].count('A')
        T_count = fa_dict[i].count('T')
        G_count = fa_dict[i].count('G')
        C_count = fa_dict[i].count('C')
        A_ratio = A_count/(A_count+T_count+C_count+G_count)
        T_ratio = T_count/(A_count+T_count+C_count+G_count)
        G_ratio = G_count/(A_count+T_count+C_count+G_count)
        C_ratio = C_count/(A_count+T_count+C_count+G_count)
        print(i,'{:.2f}'.format(A_ratio),'{:.2f}'.format(T_ratio),'{:.2f}'.format(C_ratio),'{:.2f}'.format(G_ratio),sep='\t')

if __name__ == '__main__':
    if len(sys.argv) != 2:
        Font.Scripts_tip('Statistical ATGC content in the sequence')
        print('Usage:')
        print('\tpython {0} in.fasta > stat.txt'.format(sys.argv[0]))
    else:
        Fadict = fa2dict.Fa2dict(sys.argv[1])
        count_base_number(Fadict)