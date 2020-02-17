# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

import os


os.chdir(r'D:\Result')

SpeciesID = []


# get every species id information
def get_SpeciesID(speciesID):
    with open(speciesID) as fh:
        for line in fh:
            SpeciesID.append(line.strip())

get_SpeciesID('species_id.txt')
species_list = SpeciesID
print(species_list)


with open('groups.txt') as f:
    for line in f:
        line = line.split(' ')
        # line = line.strip()
        #print(len(line))
        #print(line[1][0:5])
        #if not len(line)>13:
            # # print(SpeciesID)
            # if line[1][0:5] != line[2][0:5] :
            #     if not line[1][0:5] in line[2:]:
            #         print(line)
        for i in species_list:
            #print(i.count)
             if i.count < 1:
                 print(line)