# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-
'''
将orthofinder/orthomcl得到的直系同源文件统计直系同源ID，用于绘制Venn图(该脚本适用于5个物种的venn图)
使用：python families_count2venn_2.0.py C:\\Users\\dell\\Desktop\\ Orthogroups.tsv AT At.txt orange1 Cs.txt
     Dil Longan.txt GSVIVT Vv.txt EVM Xs.txt
此脚本是families_count2venn.py的升级版本，可直接对关键字进行设置
'''
import os
import re
import sys


# os.chdir(r'D:\Result\Transit')
os.chdir(sys.argv[1])

def famliy_count2venn(key1,key2,key3,key4,key5,key6):

    Specie1 = key1
    Specie2 = key2
    Specie3 = key3
    Specie4 = key4
    Specie5 = key5
    Specie6 = key6


    Specie1_count = open(sys.argv[4],'w')
    Specie2_count = open(sys.argv[6],'w')
    Specie3_count = open(sys.argv[8],'w')
    Specie4_count = open(sys.argv[10],'w')
    Specie5_count = open(sys.argv[12],'w')
    Specie6_count = open(sys.argv[14],'w')

    with open(sys.argv[2]) as f :
        for line in f:
            line = line.strip()
            if Specie1 in line:
                Specie1_ID = re.findall('(OG[0-9]+)',line)
                for i in Specie1_ID:
                    Specie1_count.write(i+'\n')

            if Specie2 in line :
                Specie2_ID = re.findall('(OG[0-9]+)',line)
                for i in Specie2_ID:
                    Specie2_count.write(i+'\n')

            if Specie3 in line :
                Specie3_ID = re.findall('(OG[0-9]+)',line)
                for i in Specie3_ID:
                    Specie3_count.write(i+'\n')

            if Specie4 in line:
                Specie4_ID = re.findall('(OG[0-9]+)',line)
                for i in Specie4_ID:
                    Specie4_count.write(i+'\n')

            if Specie5 in line:
                Specie5_ID = re.findall('(OG[0-9]+)',line)
                for i in Specie5_ID:
                    Specie5_count.write(i + '\n')

            if Specie6 in line:
                Specie6_ID = re.findall('(OG[0-9]+)',line)
                for i in Specie6_ID:
                    Specie6_count.write(i + '\n')

    Specie1_count.close()
    Specie2_count.close()
    Specie3_count.close()
    Specie4_count.close()
    Specie5_count.close()
    Specie6_count.close()
# 可不用 return  Specie1_count,Specie2_count,Specie3_count,Specie4_count,Specie5_count
#可用1：
#famliy_count2venn('AT','XSsp','LOC_Os','Sspon','Sobic')
#可用2：
if __name__ == '__main__':
    (famliy_count2venn(sys.argv[3],sys.argv[5],sys.argv[7],sys.argv[9],sys.argv[11],sys.argv[13]))