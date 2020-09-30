# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

#Output the fa file as a line of header and line of seq
import os
import sys

fa_dict = {}
with open(sys.argv[1]) as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            fa_key = line[1:20]
            fa_dict[fa_key] = ''
        else:
            line = line.upper()
            fa_dict[fa_key] += line.replace('\n','')
for i in fa_dict.keys():
    A_count = fa_dict[i].count('A')
    T_count = fa_dict[i].count('T')
    G_count = fa_dict[i].count('G')
    C_count = fa_dict[i].count('C')
    A_ratio = A_count/(A_count+T_count+C_count+G_count)
    T_ratio = T_count/(A_count+T_count+C_count+G_count)
    G_ratio = G_count/(A_count+T_count+C_count+G_count)
    C_ratio = C_count/(A_count+T_count+C_count+G_count)
    #print(i,'A ration is {:.2f}'.format(A_ratio),'T ration is {:.2f}'.format(T_ratio),'C ration is {:.2f}'.format(C_ratio),'G ration is {:.2f}'.format(G_ratio))
    print(i,'{:.2f}'.format(A_ratio),'{:.2f}'.format(T_ratio),'{:.2f}'.format(C_ratio),'{:.2f}'.format(G_ratio))
#step2
#for i in fa_dict.keys():
#    print(len(fa_dict[i]))   
