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
import richlog as rl


def read_geneABbed(NpXbed,AP85bed):
    header=['seqnames','start','end','ID','AB']
    Ndf = pd.read_table(NpXbed,sep='\t',names=header)
    Adf = pd.read_table(AP85bed,sep='\t',names=header)
    return Ndf,Adf


def read_anchors(anchors):
    header=['ID1','ID2','score']
    df = pd.read_table(anchors,sep='\t',comment='#',names=header)
    newID1 = df['ID1'].str.split('.t',expand=True).iloc[:,0]
    if newID1.shape[0] == df.iloc[:,[1,2]].shape[0]:
        rl.debug_out('The number of IDs after modification is the same as before modification')
        n_df = pd.concat([df.iloc[:,[1,2]],newID1],axis=1)
        n_df.rename(columns={0:'ID1'},inplace=True)
        ndf = n_df.iloc[:,[-1,0,1]]
    else:
        rl.error_out('The number of IDs before and after modification is inconsistent')
        sys.exit()
    
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
    A2B_NID = [];B2A_NID = [];Conserved_NID = []
    A2B_AID = [];B2A_AID = [];Conserved_AID = []
    for index,row in ndf.iterrows():
        if row['Nab'] == 'A' and row['Aab'] == 'B':
           A2B += 1
           A2B_NID.append(row['ID1'])
           A2B_AID.append(row['ID2'])
        elif row['Nab'] == 'B' and row['Aab'] == 'A':
            B2A += 1
            B2A_NID.append(row['ID1'])
            B2A_AID.append(row['ID2'])
        else:
            Conserved += 1
            Conserved_NID.append(row['ID1'])
            Conserved_AID.append(row['ID2'])
    return A2B,B2A,Conserved,ndf.shape[0],A2B_NID,A2B_AID,B2A_NID,B2A_AID,Conserved_NID,Conserved_AID


def output(A2B_Nid,A2B_Aid,B2A_Nid,B2A_Aid,Con_Nid,Con_Aid):
    with open('NpX2AP85.A2B.NpXgeneID.txt','w') as w1:
        for i1 in A2B_Nid:
            w1.write(i1+'\n')
    with open('NpX2AP85.A2B.AP85geneID.txt','w') as w2:
        for i2 in A2B_Aid:
            w2.write(i2+'\n')


    with open('NpX2AP85.B2A.NpXgeneID.txt','w') as w3:
        for j1 in B2A_Nid:
            w3.write(j1+'\n')
    with open('NpX2AP85.B2A.AP85geneID.txt','w') as w4:
        for j2 in B2A_Aid:
            w4.write(j2+'\n')


    with open('Conserved.NpXgeneID.txt','w') as w5:
        for k1 in Con_Nid:
            w5.write(k1+'\n')
    with open('Conserved.AP85geneID.txt','w') as w6:
        for k2 in Con_Aid:
            w6.write(k2+'\n')



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
        rl.info_out('Read the gene AB compartment bed file')
        Ndf,Adf = read_geneABbed(Nbed,Abed)
        rl.info_out('Read the anchors file')
        ndf = read_anchors(anch)
        rl.info_out('Store gene AB compartment to dictionary')
        Ngene2AB,Agene2AB = storeABinD(Ndf,Adf)
        rl.info_out('Match the AB compartment status of query gene')
        nndf = matchanchor2AB(ndf,Ngene2AB,Agene2AB)
        rl.info_out('Output the gene IDs of each species with different AB status')
        A2B,B2A,Conserved,total_otho,A2B_NID,A2B_AID,B2A_NID,B2A_AID,Conserved_NID,Conserved_AID = stat_ABswitch(nndf)
        output(A2B_NID,A2B_AID,B2A_NID,B2A_AID,Conserved_NID,Conserved_AID)
        console.print('A2B:{} , Ratio:{}'.format(A2B,round(A2B/total_otho,2)))
        console.print('B2A:{} , Ratio:{}'.format(B2A,round(B2A/total_otho,2)))
        console.print('Conserved:{} , Ratio:{}'.format(Conserved,round(Conserved/total_otho,2)))
        console.print('Total otho:{}'.format(total_otho))
    else:
        sys.exit(p.print_help())
