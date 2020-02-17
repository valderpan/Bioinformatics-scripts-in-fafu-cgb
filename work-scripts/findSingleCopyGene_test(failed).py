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

with open('groups.txt') as f:  # 共112行
    for line in f:
        line = line.strip().split()
        #print(line)
        # key = 0
        for i in species_list:
            i_count = len(re.findall(r'i',str(line)))
            if not i_count >2:
                print(line)