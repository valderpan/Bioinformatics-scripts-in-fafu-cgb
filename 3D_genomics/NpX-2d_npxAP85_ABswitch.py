#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/8/13

'''
%prog <NpXbed> <AP85bed> <anchors>

Statistics of ABswitch according to bed file (as produced by NpX-2d_npxap85_ABclassify.py) and anchors file
'''

import sys
import pandas as pd



def read_geneABbed(NpXbed,AP85bed):
    header=['seqnames','start','end','ID','AB']
    Ndf = pd.read_table(NpXbed,sep='\t',names=header)
    Adf = pd.read_table(AP85bed,sep='\t',names=header)
    return Ndf,Adf


def read_anchors(anchors):
    header=['ID1','ID2','score']
    df = pd.read_table(anchors,sep='\t',comment='#',names=header)
    newID1 = df['ID1'].str.split('.t',expand=True).iloc[:,0]
    n_df = pd.merge(df.iloc[:,[1,2]],newID1,on=df.index)
    print(n_df)
    n_df.rename(columns={0:'ID1'},inplace=True)
    ndf = n_df.iloc[:,[-1,1,2]]
    return ndf


def storeABinD(Ndf,Adf):
    Ngene2AB = Ndf.set_index('ID').to_dict()['AB']
    Agene2AB = Adf.set_index('ID').to_dict()['AB']
    return Ngene2AB,Agene2AB


def matchanchor2AB(ndf,Ngene2AB,Agene2AB):
    ndf['Nab'] = ndf['ID1'].map(Ngene2AB)
    ndf['Aab'] = ndf['ID2'].map(Agene2AB)
    return ndf


def stat_ABswitch(ndf):
    A2B = 0;B2A = 0;Conserved=0
    A2B_ID = [];B2A_ID = [];Conserved_ID = []
    for index,row in ndf.iterrows():
        if row['Nab'] == 'A' and row['Aab'] == 'B':
           A2B += 1
           A2B_ID.append(row['ID1'])
        elif row['Nab'] == 'B' and row['Aab'] == 'A':
            B2A += 1
            B2A_ID.append(row['ID2'])
        else:
            Conserved += 1
            Conserved_ID.append(row['ID2'])
    return A2B,B2A,Conserved,ndf.shape[0]


if __name__ == '__main__':
    from optparse import OptionParser
    from rich.console import Console
    from rich.traceback import install
    console = Console()
    p = OptionParser(__doc__)
    opts, args = p.parse_args()
    install()
    if len(args) == 3:
        Nbed = args[0]
        Abed = args[1]
        anch = args[2]
        Ndf,Adf = read_geneABbed(Nbed,Abed)
        ndf = read_anchors(anch)
        Ngene2AB,Agene2AB = storeABinD(Ndf,Adf)
        nndf = matchanchor2AB(ndf,Ngene2AB,Agene2AB)
        A2B,B2A,Conserved,total_otho = stat_ABswitch(nndf)
        console.print('A2B:{} , Ratio:{}'.format(A2B,round(A2B/total_otho,2)))
        console.print('B2A:{} , Ratio:{}'.format(B2A,round(B2A/total_otho,2)))
        console.print('Conserved:{} , Ratio:{}'.format(Conserved,round(Conserved/total_otho,2)))
        console.print('Total otho:{}'.format(total_otho))
    else:
        sys.exit(p.print_help())