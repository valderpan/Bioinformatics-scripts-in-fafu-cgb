# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-
import os
import re
os.chdir(r'D:\databases\Athaliana\Araport11\annotation')
with open('Athaliana_447_Araport11.protein_primaryTranscriptOnly.fa') as f :
    for line in f :
        line = line.strip()
        if line.startswith('>'):
             line = line[0:95]
             print(line)