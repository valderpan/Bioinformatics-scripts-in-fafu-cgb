# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2020/12/21

import sys
sys.path.append('../')
import pandas as pd
from path import Path
import Fontcolor as Font

pd.set_option('display.max_columns',8)

def groupby_and_concat_files(dir_path,outputfile=0):
    '''
    将plantCARE的结果进行处理并合并
    :param dir_path:
    :param outputfile:
    :return:
    '''
    all_df_list = []
    dirs = Path(dir_path).dirs()
    for dir in dirs:
        files = Path(dir).files()
        for file in files:
            if file.endswith('.tab'):
                try:
                    df = pd.read_table(file,names=['GeneID','Site Name','Sequence','Position','Length','Strand','Organism','Function'])
                    df = df.dropna(how='any')
                    df = df[~df['Function'].str.contains('\?')]
                    df = df.groupby(['GeneID', 'Site Name', 'Sequence', 'Strand', 'Length', 'Organism', 'Function'])[
                        'Position'].apply(list).reset_index()
                    all_df_list.append(df)
                except:
                    Font.Error_output('File {} read failed'.format(file))

    result = pd.concat(all_df_list,ignore_index=True)
    return result

def prepare_for_CRE_classification(CRE_df,output=0):
    '''
    将调控元件进行分类并将最终结果输出
    :param CRE_df:
    :return:
    '''
    function_dict ={'cis-acting element involved in salicylic acid responsiveness':'Salicylic acid responsiveness element',
                    '60K protein binding site':'Protein binding site',
                    'auxin-responsive element':'Auxin responsive element',
                    'binding site of AT-rich DNA binding protein (ATBP-1)':'AT-rich DNA binding protein binding site element',
                    'cis-acting element conferring high transcription levels':'High transcription levels element',
                    'cis-acting element involved in cell cycle regulation':'Cell cycle regulation element',
                    'cis-acting element involved in defense and stress responsiveness':'Defense and stress responsiveness element',
                    'cis-acting element involved in dehydration, low-temp, salt stresses':'Dehydration,low-temp,salt stresses element',
                    'cis-acting element involved in gibberellin-responsiveness':'Gibberellin-responsiveness element',
                    'cis-acting element involved in heat stress responsiveness':'Heat stress responsiveness element',
                    'cis-acting element involved in light responsiveness':'Light responsiveness element',
                    'cis-acting element involved in low-temperature responsiveness':'Low-temperature responsiveness element',
                    'cis-acting element involved in phytochrome down-regulation  expression':'Phytochrome down-regulation expression element',
                    'cis-acting element involved in the abscisic acid responsiveness':'Abscisic acid responsiveness element',
                    'cis-acting regulatory element':'Cis-acting regulatory element', #TODO
                    'cis-acting regulatory element essential for the anaerobic induction':'Anaerobic induction response element',
                    'cis-acting regulatory element involved in auxin responsiveness':'Auxin response element',
                    'cis-acting regulatory element involved in circadian control':'Circadian control response element',
                    'cis-acting regulatory element involved in endosperm regulation':'Endosperm regulation element',
                    'cis-acting regulatory element involved in light responsiveness':'Light response element',
                    'cis-acting regulatory element involved in seed-specific regulation':'Seed-specific regulation element',
                    'cis-acting regulatory element involved in the MeJA-responsiveness':'MeJA response element',
                    'cis-acting regulatory element involved in zein metabolism regulation':'Zein metabolism regulation element',
                    'cis-acting regulatory element related to meristem expression':'Zein metabolism regulation element',
                    'cis-acting regulatory element related to meristem specific activation':'Meristem specific activation element',
                    'cis-acting regulatory element root specific':'Root specific element',
                    'cis-regulatory element involved in endosperm expression':'Endosperm expression element',
                    'common cis-acting element in promoter and enhancer regions':'Common element',
                    'core promoter element around -30 of transcription start':'Core promoter element',
                    'element for maximal elicitor-mediated activation (2copies)':'Inducible promoter element',
                    'element involved in differentiation of the palisade mesophyll cells':'Mesophyll cells differentiation element',
                    'enhancer-like element involved in anoxic specific inducibility':'Anoxic inducibility enhancer-like element',
                    'gibberellin-responsive element':'Gibberellin response element',
                    'gibberellin-responsive element and part of a light responsive element':'Gibberellin,light response element',
                    'GT-1 factor binding site':'GT-1 Factor binding site',
                    'HBP factor binding site':'HBP Factor binding site',
                    'involved in activation of zein gene endosperm development':'Endosperm development element',
                    'involved in endosperm-specific negative expression':'Endosperm development element',
                    'involved in response to auxin, salicylic acid and methyl jasmonate':'Auxin,salicylic acid and methyl jasmonate response element',
                    'light responsive element':'Light response element',
                    'MeJA-responsive element':'MeJA response element',
                    'MYB binding site involved in drought-inducibility':'Drought-inducibility MYB binding site',
                    'MYB binding site involved in flavonoid biosynthetic genes regulation':'Flavonoid biosynthetic genes regulation MYB binding site',
                    'MYB binding site involved in light responsiveness':'Light response MYB binding site',
                    'MYBHv1 binding site':'MYBHv1 Binding site',
                    'part of a conserved DNA module array (CMA3)':'DNA module array conserved element',
                    'part of a conserved DNA module involved in light responsiveness':'Light response element',
                    'part of a light response element':'Light response element',
                    'part of a light responsive element':'Light response element',
                    'part of a light responsive module':'Light response element',
                    'part of a module for light response':'Light response element',
                    'part of an auxin-responsive element':'Auxin response element',
                    'part of gapA in (gapA-CMA1) involved with light responsiveness':'Light response element',
                    'protein binding site':'Protein binding site',
                    'responsive element to auxin-free medium':'Auxin response element',
                    'SEF1 factor binding site':'SEF1 Factor binding site',
                    'sequence conserved in alpha-amylase promoters':'Alpha-amylase promoters conserverd sequence element',
                    'transactivator in the cell-cycle dependent transcription':'Cell cycle regulation element',
                    'wound-responsive element':'Wound response element'
    }
    for key in function_dict:
        CRE_df.loc[CRE_df['Function'] == key,'Category']=function_dict[key]
    # print(CRE_df)
    if output:
        CRE_df.to_csv(output,sep='\t',header=True,index=False)
    return CRE_df

def unmatched_check(CREDF):
    '''
    将未进行归类(Category列为空值)的行(元件)进行输出
    :param CREDF:
    :return:
    '''
    unmatched_row = []
    for index,row in CREDF.iterrows():
        # if len(list(CREDF.loc[index].value_counts())) != 9 :
        #     unmatched_row.append(CREDF.loc[[index],:])
        if not row['Category']:
            unmatched_row.append(CREDF.loc[[index], :])
    try:
        unmatched = pd.concat(unmatched_row,ignore_index=True)
        unmatched.to_csv('03-out.unmatched.tab',sep='\t',index=False,header=True)
    except ValueError:
        Font.Log_output('All element are fully matched and have been classified')

# CRE_df = groupby_and_concat_files(r'D:\chrom download\SbSSOsZm','sctipstest.tab')
# CREDF = prepare_for_CRE_classification(CRE_df)
# unmatched_check(CREDF)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        Font.Scripts_tip('Clean,merge and classify the output from the plantcare website')
        print('Usage:')
        print('\tpython {0} <PlantcareDir.path> <out.tab>'.format(sys.argv[0]))
    else:
        CRE_df = groupby_and_concat_files(sys.argv[1],sys.argv[2])
        CREDF = prepare_for_CRE_classification(CRE_df)
        unmatched_check(CREDF)