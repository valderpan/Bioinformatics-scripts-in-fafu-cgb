# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/9/26


import sys
import pandas as pd
import math
from intervaltree import IntervalTree,Interval


def cal_Ks_valuecounts(Ks_file,limit):
    df = pd.read_table(Ks_file,sep='\t')
    ksdf = df.iloc[:,[0,3]]
    ksdf = ksdf[(ksdf['dS-ng'] >= 0) & (ksdf['dS-ng'] < limit)]
    ksdf.iloc[:,1] = round(ksdf.iloc[:,1],3)
    ks_count_dict = ksdf.iloc[:, 1].value_counts().to_dict()
    return ks_count_dict

  
def split_windows(L,S):
    limit = float(L)
    step = float(S)
    tree_D = {}
    tree = IntervalTree()
    windows = math.ceil(limit / step)
    win_s = 0
    for i in range(windows):
        if i == 0:
            win_e = win_s + step
            tree.addi(win_s, win_e, 0)
        else:
            if i < windows - 1:
                win_s += step
                if win_s + step < limit:
                    win_e = win_s + step
                else:
                    win_e = limit
                tree.addi(win_s, win_e, 0)
            else:
                win_s += step
                win_e = limit
                tree.addi(win_s, win_e, 0)
    return tree


def cal_Ks_density(ks_count_dict,tree,limit):
    for key in ks_count_dict.keys():
        inter_L = sorted(tree[key])
        if len(inter_L) == 1:
            inter = inter_L[0]
            new_inter = Interval(inter.begin, inter.end, inter.data + ks_count_dict[key])
            tree.remove(inter)
            tree.add(new_inter)
        else:
            if key == limit:
                inter = sorted(tree)[-1]
                new_inter = Interval(inter.begin,inter.end,inter.data +ks_count_dict[key])
                tree.remove(inter)
                tree.add(new_inter)
            else:
                print('This value lies in two intervals!!!')  
    return tree


def convert_result(newtree,ks_name):
    Dict = {}
    for inter in sorted(newtree):
        Dict[round(inter.begin,2)] = inter.data
    df = pd.DataFrame([Dict]).T
    df.rename(columns={0:'{}_counts'.format(ks_name)},inplace=True)
    df[ks_name] = df['{}_counts'.format(ks_name)]/sum(df['{}_counts'.format(ks_name)])*100
    output = df.sort_index().reset_index()
    return output


if __name__ == '__main__':
    output = reads_mapping_calcultate(sys.argv[1],sys.argv[2],5.0)
    output.to_excel('result1.xlsx',header=False,index=False)
    ks_count_dict = cal_Ks_valuecounts(sys.argv[1],5.0)
    tree = split_windows(5.0,0.01)
    newtree = cal_Ks_density(ks_count_dict,tree,5.0)
    o = convert_result(newtree,sys.argv[2])
    o.to_excel('result2.xlsx',header=False,index=False)