# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/9/26


import sys
import math
import argparse
import pandas as pd
from intervaltree import IntervalTree,Interval


__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20211015'
__version__ = 'v2.0'


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
    IDs = {}
    for key in ks_count_dict.keys():
        inter_L = sorted(tree[key])
        if len(inter_L) == 1:
            inter = inter_L[0]
            new_inter = Interval(inter.begin, inter.end, inter.data + ks_count_dict[key])
            if not inter.begin in IDs:
                IDs[inter.begin] = []
                IDs[inter.begin].append(key)
            else:
                IDs[inter.begin].append(key)
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
    return tree,IDs


def convert_result(newtree, ks_name):
    for inter in sorted(newtree):
        #Dict[round(inter.begin, 2)] = inter.data
        Dict[round(inter.end, 2)] = inter.data
    df = pd.DataFrame([Dict]).T
    df.rename(columns={0: '{}_counts'.format(ks_name)}, inplace=True)
    df[ks_name] = df['{}_counts'.format(ks_name)] / sum(df['{}_counts'.format(ks_name)]) * 100
    output = df.sort_index().reset_index()
    mutate  = {'index':0,'{}_counts'.format(ks_name):0,'{}'.format(ks_name):0}
    mutate_df = pd.DataFrame([mutate])
    output = pd.concat([mutate_df,output])
    return output


def main(args):
    from rich.traceback import install
    install()
    KaKs_res = args.kaks
    limit = args.limit
    window = args.window
    KaKs_name = args.name
    ks_count_D = cal_Ks_valuecounts(KaKs_res,limit)
    Tree = split_windows(limit,window)
    newtree,IDs_D = cal_Ks_density(ks_count_D,Tree,limit)
    output = convert_result(newtree,KaKs_name)
    print(output)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''calculate Ks distribution''',
        usage="python {} -k SsSo.KaKs.result -n SsSo -l limit -w window > output.bed".format(sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,
                                                                            __mail__, __date__, __version__))
    parser.add_argument('-k', '--kaks', required=True, help='Input KaKs result(.KaKs.result)')
    parser.add_argument('-l', '--limit', required=True, type=float,help='Input Ks limit[float]')
    parser.add_argument('-w', '--window', required=True,type=float, help='Input window size[float]')
    parser.add_argument('-n', '--name', required=True,type=float, help='Input KaKs_file name')
    args = parser.parse_args()
    main(args)
