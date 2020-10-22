#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan



import sys
import pandas as pd
import re
import numpy as np
import argparse

__author__ = 'Haoran Pan'
__mail__ = 'haoran_pan@qq.com'
__date__ = '20201022'
__version__ = 'v1.0'

def sort_bed(bed_file):
    bed_df = pd.read_table(bed_file,sep='\t',names=['chrid','start','end','geneid','tmp','strand'])
    chr_id = bed_df['chrid'].unique()
    concat_list = []
    for i in chr_id:
        bed_chr_df = bed_df[bed_df['chrid'].str.contains(i)]
        bed_chr_df = bed_chr_df.sort_values('start')
        concat_list.append(bed_chr_df)
    result = pd.concat(concat_list,ignore_index=True)
    result = result[~result['chrid'].str.contains('scaffold', 'contig')]
    bed_file_name = re.findall('([0-9A-Za-z]+)\.bed',bed_file)[0]
    result.to_csv(bed_file_name+'_sorted.bed',sep='\t',header=None,index=False)
    return result


def block_complete(bed_file,anchors_file):
    bed_df = sort_bed(bed_file)
    anchors_df = pd.read_table(anchors_file,sep='\t',names=['speciesA','speciesB','score'])

    # ----------按照'###'将blocks进行拆分----------
    split_row_list = []
    for index,row in anchors_df.iterrows():
        if anchors_df.loc[index]['speciesA'] == '###':
            split_row_list.append(index)
    concat_row_list = []
    for i in range(0,len(split_row_list)):
        if i != len(split_row_list) -1 :
            blocks = anchors_df.iloc[split_row_list[i]:split_row_list[i+1],:]
            concat_row_list.append(blocks)
        else:
            blocks = anchors_df.iloc[split_row_list[i]:anchors_df.shape[0],:]
            concat_row_list.append(blocks)
    #------------------------------------------
    inter_blocks_list = []
    for i in concat_row_list:
        first_geneid = i.iloc[1,:]['speciesA'];second_geneid = i.iloc[2,:]['speciesA'];third_geneid = i.iloc[3,:]['speciesA']
        last_geneid = i.iloc[-1, :]['speciesA'];second_last_geneid = i.iloc[-2, :]['speciesA'];third_last_geneid = i.iloc[-3,:]['speciesA']
        
        #----------将每个block的首尾基因之间的基因在bed文件中找到，对block进行补全----------
        first_row = bed_df[bed_df['geneid'] == first_geneid].index.tolist()
        if len(first_row)>0:
            first_row_index = bed_df[bed_df['geneid'] == first_geneid].index.tolist()[0]
        else:
            second_row = bed_df[bed_df['geneid'] == second_geneid].index.tolist()
            if len(second_row) >0:
                first_row_index = bed_df[bed_df['geneid'] == second_geneid].index.tolist()[0]
            else:
                third_row = bed_df[bed_df['geneid'] == third_geneid].index.tolist()
                if len(third_row) > 0:
                    first_row_index = bed_df[bed_df['geneid'] == third_geneid].index.tolist()[0]
                else:
                    continue
        last_row = bed_df[bed_df['geneid'] == last_geneid].index.tolist()
        if len(last_row)>0:
            last_row_index = bed_df[bed_df['geneid'] == last_geneid].index.tolist()[0]
        else:
            second_last_row = bed_df[bed_df['geneid'] == second_last_geneid].index.tolist()
            if len(second_last_row) > 0:
                last_row_index = bed_df[bed_df['geneid'] == second_last_geneid].index.tolist()[0]
            else:
                third_last_row = bed_df[bed_df['geneid'] == third_last_geneid].index.tolist()
                if len(third_last_row)>0:
                    last_row_index = bed_df[bed_df['geneid'] == third_last_geneid].index.tolist()[0]
                else:
                    continue
        #--------------------------------------------------------------------------

        blocks_rowindex_list = [num for num in range(first_row_index,last_row_index+1)]

        #----------存储每个区块的geneid并根据geneid在anchors文件中进行匹配----------
        blocks_geneid_list = []
        for index,row in bed_df.iterrows():
            if index in blocks_rowindex_list:
                blocks_geneid_list.append(bed_df.loc[index]['geneid'])

        intra_blocks_concat = []
        for gene in blocks_geneid_list:
            if gene in i.iloc[:,0].tolist():
                intra_blocks_concat.append(i[i.loc[:,'speciesA'] == gene])
            else:
                tmp = {'speciesA':gene,'speciesB':np.nan,'score':np.nan}
                intra_blocks_concat.append(pd.DataFrame(tmp, index=[0]))

        intra_blocks_result = pd.concat(intra_blocks_concat,ignore_index=True)
        intra_blocks_result.loc[-1] = ['###',np.nan,np.nan]
        intra_blocks_result = intra_blocks_result.sort_index()
        inter_blocks_list.append(intra_blocks_result)
        #-------------------------------------------------------------------
    inter_blocks_result = pd.concat(inter_blocks_list,ignore_index=True)
    #print(inter_blocks_result)
    #inter_blocks_result.to_csv('01_blocks_complete.tab',sep='\t',header=None,index=False)
    return inter_blocks_result


def block_match(ANC_file,block_complete_file):
    AMK = pd.read_table(ANC_file,sep='\t')
    AMK = AMK.dropna()
    AMK_Dict = AMK.set_index('Gene').to_dict()['ancestor_chromosome']
    # blocks_supplement = pd.read_table(block_complete_file,sep='\t',names=['speciesA','speciesB','score'])
    blocks_supplement = block_complete_file
    # ----------按照'###'将blocks进行拆分----------
    split_row_list = []
    for index, row in blocks_supplement.iterrows():
        if blocks_supplement.loc[index]['speciesA'] == '###':
            split_row_list.append(index)
    concat_row_list = []
    for i in range(0, len(split_row_list)):
        if i != len(split_row_list) - 1:
            blocks = blocks_supplement.iloc[split_row_list[i]:split_row_list[i + 1], :]
            concat_row_list.append(blocks)
        else:
            blocks = blocks_supplement.iloc[split_row_list[i]:blocks_supplement.shape[0], :]
            concat_row_list.append(blocks)
    # ------------------------------------------
    inter_blocks_list = []
    for block in concat_row_list:
        block['Chromosome'] = block['speciesB'].map(AMK_Dict)
        if len(block['Chromosome'].value_counts().tolist())>0:
            inter_blocks_list.append(block.iloc[:,[0,3]])
        else:
            inter_blocks_list.append(block.iloc[[0],:].iloc[:,[0,3]])
    inter_blocks_result = pd.concat(inter_blocks_list)
    # print(inter_blocks_result)
    inter_blocks_result.to_csv('block_match.tab',sep='\t',header=None,index=False)
    # return inter_blocks_result


def bed_map(block_match_file,output_file):
    block_df = pd.read_table(block_match_file,sep='\t',names=['speciesA','Chromosome'])
    # block_df = block_match_file

    split_row_list = []
    for index, row in block_df.iterrows():
        if block_df.loc[index]['speciesA'] == '###':
            split_row_list.append(index)
    concat_row_list = []
    for i in range(0, len(split_row_list)):
        if i != len(split_row_list) - 1:
            blocks = block_df.iloc[split_row_list[i]:split_row_list[i + 1], :]
            concat_row_list.append(blocks)
        else:
            blocks = block_df.iloc[split_row_list[i]:block_df.shape[0], :]
            concat_row_list.append(blocks)

    inter_blocks_list = []
    for block in concat_row_list:
        if block.shape[0] > 1:
            block['Chromosome'] = block['Chromosome'].value_counts().index[0]
            block = block[~block['speciesA'].str.contains('###')]
            inter_blocks_list.append(block)
    inter_blocks_result = pd.concat(inter_blocks_list)
    # print(inter_blocks_result)
    inter_blocks_result.to_csv(output_file,sep='\t',header=None,index=False)


def main(args):
    if args.bedfile and args.anchors:
        bed_file = args.bedfile
        anchors_file = args.anchors
        ANC_file = args.Ancestralchr
        block_complete_file = block_complete(bed_file, anchors_file)
        block_match(ANC_file, block_complete_file)
    elif args.blockmatch and args.output:
        block_match_file = args.blockmatch
        output = args.output
        bed_map(block_match_file, output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''test script''',
        epilog= 'author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
    parser.add_argument('-b', '--bedfile', help='Input jcvi syntenic analysis file. eg:speciesA.bed')
    parser.add_argument('-a', '--anchors', help='Input the anchors file of the jcvi syntenic analysis result. eg:speciesA.speciesB.anchors')
    parser.add_argument('-A', '--Ancestralchr', help='Input the match file of the ancestor chromosome of monocot/dicot. [monocot:AMK.txt]/[Dicots:AEK.txt]')

    parser.add_argument('-m', '--blockmatch', help='Intput the block_match file generated in the step1')
    parser.add_argument('-o', '--output', help='Output the final map file .Suggested naming: Anc.map')
    args = parser.parse_args()
    main(args)
