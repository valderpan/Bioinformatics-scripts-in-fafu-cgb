# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-


## 将单拷贝基因的基因ID找出来，方便下一步去每个物种的蛋白序列中根据ID找序列
import os
import re
#
#
os.chdir(r'D:\Result')


Ptri_list=[]
Tcac_list=[]
Vvin_list=[]
Zmay_list=[]
long_list=[]
Macu_list=[]
Osat_list=[]
Ccle_list=[]
Cpap_list=[]
Atri_list=[]
Atha_list=[]
Lchi_list=[]
ramb_list=[]

with open('SingleCopyGene.txt') as f:
    for line in f:
        line = line.strip()

        line =line.split(", ")

        for i in line:
            if i[1:5] == "Ptri":
                Ptri_list.append(i)
            if i[1:5] == 'Tcac':
                Tcac_list.append(i)
            if i[1:5] == 'Vvin':
                Vvin_list.append(i)
            if i[1:5] == 'Zmay':
                Zmay_list.append(i)
            if i[1:5] == 'long':
                long_list.append(i)
            if i[1:5] == 'Macu':
                Macu_list.append(i)
            if i[1:5] == 'Osat':
                Osat_list.append(i)
            if i[1:5] == 'Ccle':
                Ccle_list.append(i)
            if i[1:5] == 'Cpap':
                Cpap_list.append(i)
            if i[1:5] == 'Atri':
                Atri_list.append(i)
            if i[1:5] == 'Atha':
                Atha_list.append(i)
            if i[1:5] == 'Lchi':
                Lchi_list.append(i)
            if i[1:5] == 'ramb':
                ramb_list.append(i)

    print(Ptri_list)
    print(Tcac_list)
    print(Vvin_list)
    print(long_list)
    print(Macu_list)
    print(Osat_list)
    print(Ccle_list)
    print(Cpap_list)
    print(Atha_list)
    print(Atri_list)
    print(Zmay_list)
    print(Lchi_list)
    print(ramb_list)
'''
上一步输出结果为：
["'Athaliana|AT3G12650.1']", "'Athaliana|AT2G26070.2']"...
下一步将结果中的所有标点去掉
'''

## 运行下面代码最好将上边不必要部分#掉

with open('scp_for_each_species_1.txt') as tmp_f:
    for line in tmp_f:
        line = line.rstrip()
        line_result = re.findall(r"[^\'\[\]\"\(\)]",line)
        print(''.join(line_result))
