#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/6/22

'''
%prog <ZG_correct_round2.tigmovestat.allchr.xlsx> <movestat.concat.clean.xlsx>

将每一条染色体的移动进行合并，合并为一个总表，对其进行去冗余。
因为会出现一个contig从一个group移动到两个group中的情况。
'''
__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20210622'
__version__ = 'v1.0'

import pandas as pd
import re
from collections import Counter
import argparse
import sys
sys.path.append(r'D:\pycharm\panhr\FAFU_CGB\Genomic_exercise')
import SetLog



def need2movetig(Allchr_tigmovestat_xlsx, RefChrs,output=0):
    '''
    读取所有染色体的tigmovestat,每条染色体为一个sheet,将其合并为一个表格并输出
    :param Allchr_tigmovestat_xlsx:
    :param RefChrs: 参考物种的染色体数目
    :return:
    '''
    concat_list = []
    for num in range(0, int(RefChrs)):
        df = pd.read_excel(Allchr_tigmovestat_xlsx, sheet_name=num, header=None)
        if df.shape[0] > 0:
            df = df.dropna(how='any')
            concat_list.append(df)
        else:
            SetLog.error_out('Chr{} is empty !'.format(num + 1))
    movestat = pd.concat(concat_list, ignore_index=True)
    SetLog.info_out('Origin movestat num : {}'.format(movestat.shape[0]))

    staequalact = []  # 移动前后都是同一个group的情况
    for index, row in movestat.iterrows():
        if row[1] == row[2]:
            staequalact.append(index)
    SetLog.info_out('staequalact num : {}'.format(len(staequalact)))

    movestat = movestat[~movestat.index.isin(staequalact)]
    SetLog.info_out('removed stay contig num : {}'.format(movestat.shape[0]))
    if output:
        movestat.to_excel('movestat.concat.xlsx', header=False, index=False)
    # print(movestat)
    return movestat


def read_file(movestat,gmap,movestat_xlsx=0):
    if movestat_xlsx:
        df = pd.read_excel(movestat_xlsx, names=['tigID', 'pos1', 'pos2'])
    else:
        df = movestat
        df.rename(columns={0:'tigID',1:'pos1',2:'pos2'},inplace=True)
    # print(df)
    gmap_df = pd.read_table(gmap,comment='#',header=None)
    return df,gmap_df


def movestat_clean(df):
    '''
    找出movestat中的重复移动的contig
    :param df:
    :return:
    '''
    dup_df = df[df['tigID'].duplicated()]
    # print(dup_df)
    dup_contig = dup_df['tigID'].tolist()
    return dup_df,dup_contig


def parse_gmap(gmap_df):
    gdf = gmap_df.iloc[:,[0,2,8]]
    gdf = gdf[gdf[2] == 'gene']
    gdf_geneID = gdf[8].str.split(';',expand=True).iloc[:,0].str.split("=",expand=True).iloc[:,1]
    gdf['new'] = gdf_geneID
    gdf = gdf.iloc[:,[0,-1]]
    # print(gdf)
    tig2chr = {}
    tig2bestchr = {}

    for index,row in gdf.iterrows():
        if row[0] not in tig2chr:
            # print(row['new'])
            if 'K' not in row['new']: #会有Sobic.K009800.3.path1这类基因
                chrID  = re.findall('\.([0-9]+)G',row['new'])[0][1:]
                tig2chr[row[0]] = ['Chr{}'.format(chrID)]
        else:
            if 'K' not in row['new']:
                chrID = re.findall('\.([0-9]+)G', row['new'])[0][1:]
                tig2chr[row[0]].append('Chr{}'.format(chrID))


    for key in tig2chr.keys():
        # print(key)
        # print(tig2chr[key])
        if len(tig2chr[key]) == 1:
            tig2bestchr[key] = tig2chr[key][0]
        else:
            tmp = Counter(tig2chr[key])
            # print(tmp)
            # print(tmp.most_common(1)[0][0])
            tig2bestchr[key] = tmp.most_common(1)[0][0]
    print(len(tig2bestchr))

    return gdf,tig2bestchr


def Chr2stablegroup(refChrNum):
    Chr2groupDict = {}
    for chr in range(1,refChrNum+1):
        if chr == 1:
            Chr2groupDict['Chr0{}'.format(chr)] = ['group'+ str(x) for x in [4,7,14,31,37,47,52,65]]
        elif chr == 2:
            Chr2groupDict['Chr0{}'.format(chr)] = ['group'+ str(x) for x in [3,6,15,20,21,25,54,63]]
        elif chr == 3:
            Chr2groupDict['Chr0{}'.format(chr)] = ['group'+ str(x) for x in [1,11,12,19,22,24,59,73]]
        elif chr == 4:
            Chr2groupDict['Chr0{}'.format(chr)] = ['group'+ str(x) for x in [17,23,26,35,40,57,62,69]]
        elif chr == 5:
            Chr2groupDict['Chr0{}'.format(chr)] = ['group'+ str(x) for x in [36,39,43,46,50,68,71,77]]
        elif chr == 6:
            Chr2groupDict['Chr0{}'.format(chr)] = ['group'+ str(x) for x in [5,18,44,51,61,67,70,72]]
        elif chr == 7:
            Chr2groupDict['Chr0{}'.format(chr)] = ['group'+ str(x) for x in [9,10,28,33,34,66,79,80]]
        elif chr == 8:
            Chr2groupDict['Chr0{}'.format(chr)] = ['group'+ str(x) for x in [2,27,49,53,55,56,75,78]]
        elif chr == 9:
            Chr2groupDict['Chr0{}'.format(chr)] = ['group'+ str(x) for x in [8,13,30,32,38,60,64,74]]
        elif chr == 10:
            Chr2groupDict['Chr{}'.format(chr)] = ['group'+ str(x) for x in [16,29,41,42,45,48,58,76]]
    return Chr2groupDict


def match_contig(dup_contig,chrgroupdf,tig2bestchr,df):
    remove_Index = []
    for ctg in dup_contig:
        # print(ctg)
        # print(tig2bestchr[ctg])
        # print(chrgroupdf[tig2bestchr[ctg]])
        query_ctg = df[df['tigID']==ctg]
        # print(query_ctg)

        tmp = []
        for index,row in query_ctg.iterrows():
            if row['pos2'] not in chrgroupdf[tig2bestchr[ctg]]:
                # remove_Index.append(index)
                tmp.append(index)
        if len(tmp) == 1:
            remove_Index.extend(tmp)
        else:
            continue

        # print(remove_Index)

        # for groupID in chrgroupdf[tig2bestchr[ctg]]:
        #     if groupID not in query_ctg_pos2:
        #         remove_ID.append()
        # print(query_ctg)

        # break
    return remove_Index


def remove_dup_row(df,remove_index,result_xlsx):
    # print(df.shape[0])
    # print(len(remove_index))
    result = df[~df.index.isin(remove_index)]
    # print(result)
    result = result.drop_duplicates('tigID',keep='first')
    # result.to_excel(r'E:\Badila\synteny_analysis\Badila染色体第二次校正\矫正后\movestat_round2_cooncat.clean.xlsx',header=None,index=False)
    result.to_excel(result_xlsx,header=None,index=False)


def main(args):
    allchrmovestat = args.allstat
    gmap = args.gmap
    refchrs = args.refchrs
    result_xlsx = args.result
    concat_output = args.concatoutput
    movestat = need2movetig(allchrmovestat,refchrs,concat_output)
    df, gmap_df = read_file(movestat, gmap)
    dup_df, dup_contig = movestat_clean(df)
    gdf, tig2bestchr = parse_gmap(gmap_df)
    chrgroupdf = Chr2stablegroup(refchrs)
    remove_index = match_contig(dup_contig, chrgroupdf, tig2bestchr, df)
    remove_dup_row(df, remove_index,result_xlsx)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''The movement of each chromosome is merged, merged into a total table, and de-redundant it.\nBecause there will be a situation where a contig moves from one group to two groups.''',
        usage="python {} -a allchr.movestat.xlsx -g gmap.gff3 -r 10 -R movestat.clean.xlsx (-c 1)".format(sys.argv[0]),
        epilog= 'author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
    parser.add_argument('-a','--allstat',required=True,help='Input allchr.movestat.xlsx')
    parser.add_argument('-g','--gmap',required=True,help='Input gmap.gff3')
    parser.add_argument('-r','--refchrs',type=int,required=True,help='Input Ref species chr numbers')
    parser.add_argument('-R','--result',required=True,help='Output the result file(.xlsx)')
    parser.add_argument('-c','--concatoutput',required=False,
                help='If specified, it will output the movement table of all chromosomes congtig(movestat.concat.xlsx)')
    args = parser.parse_args()
    main(args)

# if __name__ == '__main__':
#     movestat = need2movetig(r'E:\Badila\synteny_analysis\Badila染色体第二次校正\矫正后\ZG_correct_round2.tigmovestat.allchr.xlsx',10)
#     df,gmap_df = read_file(movestat,r'E:\Badila\synteny_analysis\Badila染色体第二次校正\矫正后\gmap.gff3')
#     dup_df,dup_contig = movestat_clean(df)
#     gdf,tig2bestchr = parse_gmap(gmap_df)
#     chrgroupdf = Chr2stablegroup(10)
#     remove_index = match_contig(dup_contig,chrgroupdf,tig2bestchr,df)
#     remove_dup_row(df,remove_index)