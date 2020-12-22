#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2020/12/22

import sys
sys.path.append('../')
from intervals import IntInterval
import gain_geneCRE as geneCRE
import Fontcolor as Font

def stateCREdistribution(querygeneCRE,outputfile):
    '''
    统计每个gene上CRE在0-2000bp上20个bin上的位置
    :param querygeneCRE:
    :param outputfile:
    :return:
    '''
    data_list = [0, 100,200,300,400,500,600,700,800,900,1000,
                 1100,1200,1300,1400,1500,1600,1700,1800,1900,2000]
    list_interval = []
    # 获得两两之间的区间范围
    for i in range(len(data_list)):
        if i != len(data_list) - 1:  # 为了排除最后一个位置，索引超出矩阵维度
            data_range = IntInterval.closed(data_list[i], data_list[i + 1])  # 这里用的双闭[2, 3],两边都是闭合区间，可以用open_closed()→(],或者closed_open()→[)
            list_interval.append(data_range)
    # print(list_interval)
    output = open(outputfile,'w')
    with open(querygeneCRE) as f:
        lines = (line.strip() for line in f)
        output.write('GeneID'+'\t'+'Site_name'+'\t'+'bin1'+'\t'+'bin2'+'\t'+'bin3'+'\t'+'bin4'+'\t'+'bin5'+'\t'+'bin6'+'\t'+'bin7'+'\t'+
                     'bin8'+'\t'+'bin9'+'\t'+'bin10'+'\t'+'bin11'+'\t'+'bin12'+'\t'+'bin13'+'\t'+'bin14'+'\t'+'bin15'+'\t'+'bin16'+'\t'+'bin17'+'\t'+
                     'bin18'+'\t'+'bin19'+'\t'+'bin20'+'\n')
        for line in lines:
            Group1 = 0;Group2 = 0;Group3 = 0;Group4 = 0;Group5 = 0;Group6 = 0;Group7 = 0;Group8 = 0;Group9 = 0
            Group10 = 0;Group11 = 0;Group12 = 0;Group13 = 0;Group14 = 0;Group15 = 0;Group16 = 0;Group17 = 0;Group18 = 0
            Group19 = 0;Group20 = 0

            location_query_interval = []
            if line.startswith('GeneID'):
                continue
            else:
                line_list = line.split('\t')
                geneID = line_list[0];Sitename = line_list[1];Position = line_list[-1]
                for single_interval in list_interval:
                    for i in Position[1:-1].split(','):
                        if i in single_interval:
                            location_query_interval.append(single_interval)
            # print(location_query_interval)
            for j in location_query_interval:
                if j.lower==0 and j.upper == 100:
                    Group1 +=1
                elif j.lower==100 and j.upper == 200:
                    Group2 +=1
                elif j.lower==200 and j.upper == 300:
                    Group3 +=1
                elif j.lower==300 and j.upper == 400:
                    Group4 +=1
                elif j.lower==400 and j.upper == 500:
                    Group5 +=1
                elif j.lower==500 and j.upper == 600:
                    Group6 +=1
                elif j.lower==600 and j.upper == 700:
                    Group7 +=1
                elif j.lower==700 and j.upper == 800:
                    Group8 +=1
                elif j.lower==800 and j.upper == 900:
                    Group9 +=1
                elif j.lower==900 and j.upper == 1000:
                    Group10 +=1
                elif j.lower==1000 and j.upper == 1100:
                    Group11 +=1
                elif j.lower==1100 and j.upper == 1200:
                    Group12 +=1
                elif j.lower==1200 and j.upper == 1300:
                    Group13 +=1
                elif j.lower==1300 and j.upper == 1400:
                    Group14 +=1
                elif j.lower==1400 and j.upper == 1500:
                    Group15 +=1
                elif j.lower==1500 and j.upper == 1600:
                    Group16 +=1
                elif j.lower==1600 and j.upper == 1700:
                    Group17 +=1
                elif j.lower==1700 and j.upper == 1800:
                    Group18 +=1
                elif j.lower==1800 and j.upper == 1900:
                    Group19 +=1
                elif j.lower==1900 and j.upper == 2000:
                    Group20 +=1

            output.write(geneID+'\t'+Sitename+'\t'+str(Group1)+'\t'+str(Group2)+'\t'+str(Group3)+'\t'+str(Group4)+'\t'+str(Group5)+'\t'+
                         str(Group6)+'\t'+str(Group7)+'\t'+str(Group8)+'\t'+str(Group9)+'\t'+str(Group10)+'\t'+str(Group11)+'\t'+
                         str(Group12)+'\t'+str(Group13)+'\t'+str(Group14)+'\t'+str(Group15)+'\t'+str(Group16)+'\t'+str(Group17)+'\t'+
                         str(Group18)+'\t'+str(Group19)+'\t'+str(Group20)+'\n')

    output.close()


if __name__ == '__main__':
    if len(sys.argv) != 4:
        Font.Scripts_tip('Count the distribution of CRE on each gene')
        print('Usage:')
        print('\tpython {0} <geneID.txt> <allgenomeCRE> <CREdisOut.bed>'.format(sys.argv[0]))
    else:
        querygeneCRE = geneCRE(sys.argv[1],sys.argv[2])
        stateCREdistribution(querygeneCRE,sys.argv[3])