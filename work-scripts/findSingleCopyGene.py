# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

import os
import re
os.chdir(r'D:\Result')


SpeciesID = []

def get_SpeciesID(speciesID):
    with open(speciesID) as fh:
        for line in fh:
            SpeciesID.append(line.strip())

get_SpeciesID('species_id.txt')
species_list = SpeciesID
print(species_list)

## 打开得到的groups.txt文件，将单拷贝基因找出来

with open('groups.txt') as f:
    for line in f:
        line = line.strip().split()
        #print(line)
        # key = 0
        # for i in species_list:
        #     i_count = len(re.findall(r'i',str(line)))
        #     if not i_count >2:
        #         print(line)
        Atr_counts = len(re.findall(r'Atrichopoda', str(line)))
        Vvi_counts = len(re.findall(r'Vvinifera', str(line)))
        Ccl_counts = len(re.findall(r'Cclementina', str(line)))
        Ath_counts = len(re.findall(r'Athaliana', str(line)))
        Cpa_counts = len(re.findall(r'Cpapaya', str(line)))
        Mac_counts = len(re.findall(r'Macuminata', str(line)))
        Osa_counts = len(re.findall(r'Osativa', str(line)))
        Tca_counts = len(re.findall(r'Tcacao', str(line)))
        Lon_counts = len(re.findall(r'longyan', str(line)))
        Zma_counts = len(re.findall(r'Zmays', str(line)))
        Ptr_counts = len(re.findall(r'Ptrichocarpa',str(line)))
        Lch_counts = len(re.findall(r'Lchinese',str(line)))
        Ram_counts = len(re.findall(r'rambutan',str(line)))
        # print(Ptr_counts)
        # print(Zma_counts)
        if not Atr_counts>1:
            if not Vvi_counts > 1:
                if not Ccl_counts > 1:
                    if not Ath_counts > 1:
                        if not Cpa_counts > 1:
                            if not Mac_counts > 1:
                                if not Osa_counts > 1 :
                                    if not Tca_counts > 1 :
                                        if not Lon_counts > 1:
                                            if not Zma_counts >1 :
                                                if not Ptr_counts > 1:
                                                    if not Lch_counts >1 :
                                                        if not  Ram_counts > 1:
                                                            print(line)

