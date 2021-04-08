#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/4/8

import sys
from path import Path
from collections import OrderedDict
sys.path.append('/share/home/stu_panhaoran/scripts')
import Fontcolor as Font


def findresultfile(reslut_path):

    resultdir = [dir for dir in Path(reslut_path).dirs() if 'hic_results' in dir][0]
    stats_dir  = [dir for dir in Path(resultdir).dirs() if 'stats' in dir][0]
    sample_dirs = [dir for dir in Path(stats_dir).dirs()]
    return sample_dirs


def statsresult(dirlist):

    for sample in dirlist:
        item2num = OrderedDict()
        files = Path(sample).files()
        for file in files:
            if file.endswith('mpairstat'):
                with open(file,'r') as f1:
                    lines = (line.strip() for line in f1)
                    for line in lines:
                        if line.startswith('#'):
                            continue
                        else:
                            line_list = line.split('\t')
                            if line_list[0] == 'Total_pairs_processed':
                                item2num['Clean Paired-end Reads'] = line_list[1]
                            elif line_list[0] == 'Unmapped_pairs':
                                item2num['Unmapped Paired-end Reads'] = line_list[1]
                                item2num['Unmapped Paired-end Reads Rate (%)'] = line_list[2]
                            elif line_list[0] == 'Pairs_with_singleton':
                                item2num['Paired-end Reads with Singleton'] = line_list[1]
                                item2num['Paired-end Reads with Singleton Rate(%)'] = line_list[2]
                            elif line_list[0] == 'Multiple_pairs_alignments':
                                item2num['Multi Mapped Paired-end Reads'] = line_list[1]
                                item2num['Multi Mapped Ratio (%)'] = line_list[2]
                            elif line_list[0] == 'Unique_paired_alignments':
                                item2num['Unique Mapped Paired-end Reads'] = line_list[1]
                                item2num['Unique Mapped Ratio (%)'] = line_list[2]
            elif file.endswith('.mRSstat'):
                with open(file,'r') as f2:
                    lines = (line.strip() for line in f2)
                    for line in lines:
                        if line.startswith('#'):
                            continue
                        else:
                            line_list = line.split('\t')
                            if line_list[0] == 'Dangling_end_pairs':
                                item2num['Dangling End Paired-end Reads'] = line_list[1]
                            elif line_list[0] == 'Self_Cycle_pairs':
                                item2num['Self Circle Paired-end Reads'] = line_list[1]
                            elif line_list[0] == 'Dumped_pairs':
                                item2num['Dumped Paired-end Reads'] = line_list[1]
                            elif line_list[0] == 'Valid_interaction_pairs':
                                item2num['Interaction Paired-end Reads'] = line_list[1]
            elif file.endswith('.mergestat'):
                with open(file,'r') as f3:
                    lines = (line.strip() for line in f3)
                    for line in lines:
                        if line.startswith('#'):
                            continue
                        else:
                            line_list = line.split('\t')
                            if line_list[0] == 'valid_interaction_rmdup':
                                item2num['Lib Valid Paired-end Reads'] = line_list[1]
        print(sample+':')
        print('Statistics of mapping')
        print('\t','Clean Paired-end Reads',item2num['Clean Paired-end Reads'],sep='\t')
        print('\t','Unmapped Paired-end Reads',item2num['Unmapped Paired-end Reads'],sep='\t')
        print('\t','Unmapped Paired-end Reads Rate (%)',item2num['Unmapped Paired-end Reads Rate (%)'],sep='\t')
        print('\t','Paired-end Reads with Singleton',item2num['Paired-end Reads with Singleton'],sep='\t')
        print('\t','Paired-end Reads with Singleton Rate(%)',item2num['Paired-end Reads with Singleton Rate(%)'],sep='\t')
        print('\t','Multi Mapped Paired-end Reads',item2num['Multi Mapped Paired-end Reads'],sep='\t')
        print('\t','Multi Mapped Ratio (%)',item2num['Multi Mapped Ratio (%)'],sep='\t')
        print('\t','Unique Mapped Paired-end Reads',item2num['Unique Mapped Paired-end Reads'],sep='\t')
        print('\t','Unique Mapped Ratio (%)',item2num['Unique Mapped Ratio (%)'],sep='\t')

        print('Statistics of valid reads')
        print('\t','Unique Mapped Paired-end Reads',item2num['Unique Mapped Paired-end Reads'],sep='\t')
        print('\t','Dangling End Paired-end Reads',item2num['Dangling End Paired-end Reads'],sep='\t')
        print('\t', 'Dangling End Rate (%)',round(int(item2num['Dangling End Paired-end Reads']) / int(item2num['Unique Mapped Paired-end Reads']) *100, 3), sep='\t')
        print('\t', 'Self Circle Paired-end Reads', item2num['Self Circle Paired-end Reads'], sep='\t')
        print('\t', 'Self Circle Rate (%)',round(int(item2num['Self Circle Paired-end Reads']) / int(item2num['Unique Mapped Paired-end Reads']) *100, 3), sep='\t')
        print('\t', 'Dumped Paired-end Reads', item2num['Dumped Paired-end Reads'], sep='\t')
        print('\t', 'Dumped Rate (%)', round(int(item2num['Dumped Paired-end Reads']) / int(item2num['Unique Mapped Paired-end Reads']) *100, 3),sep='\t')
        print('\t', 'Interaction Paired-end Reads', item2num['Interaction Paired-end Reads'], sep='\t')
        print('\t', 'Interaction Rate (%)',round(int(item2num['Interaction Paired-end Reads']) / int(item2num['Unique Mapped Paired-end Reads']) *100 , 3), sep='\t')
        print('\t', 'Lib Valid Paired-end Reads', item2num['Lib Valid Paired-end Reads'], sep='\t')
        Lib_Valid_Rate = round(int(item2num['Lib Valid Paired-end Reads']) / int(item2num['Interaction Paired-end Reads']) *100, 3)
        print('\t', 'Lib Valid Rate (%)',Lib_Valid_Rate, sep='\t')

        print('\t','Lib Dup (%)',round(100-Lib_Valid_Rate,3),sep='\t')


if __name__ == '__main__':
    if __name__ == '__main__':
        if len(sys.argv) != 2:
            Font.Scripts_tip('Convert HiC-Pro output stats results into publishable level tables\n'
                             'Note: This script is available in both absolute and relative paths')
            print('Usage:')
            print('\tpython {0} hicpro_output(HiC-Pro -o parameter specifies the directory)'.format(sys.argv[0]))
        else:
            dir_list = findresultfile(sys.argv[1])
            statsresult(dir_list)