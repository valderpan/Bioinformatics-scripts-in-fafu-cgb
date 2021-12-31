#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# @Time : 2020/8/15 21:43

import pandas as pd
import numpy as np
import sys
import os
import re
from path import Path

### ----------------2020.08.30 统计NBSgene中在Nepal到AP85演化过程中compartment状态的转变-----------
# with open(r'C:\\Users\\dell\\Desktop\\Nepal_project\\20200829 NBS_gene\\1-NBS-Nipor-1612-list.txt') as f:
#     # for line in f:
#     #     print(line.strip())
#     NBS = [line.strip() for line in f ]
#     print(len(NBS))
# with open(r'C:\Users\dell\Desktop\Nepal_project\20200817 four_species_compartment_conserved_GOenrich\Himalaya.AP85.A2B_compartment.HimalayageneID.txt') as f1:
#     A2B = [line.strip() for line in f1]
#     print(len(A2B))
# with open(r'C:\Users\dell\Desktop\Nepal_project\20200817 four_species_compartment_conserved_GOenrich\Himalaya.AP85.B2A_compartment.HimalayageneID.txt') as f2:
#     B2A = [line.strip() for line in f2]
#     print(len(B2A))
# with open(r'C:\Users\dell\Desktop\Nepal_project\20200817 four_species_compartment_conserved_GOenrich\Himalaya.AP85.conserved_compartment.HimalayageneID.txt') as f3:
#     conserved = [line.strip() for line in f3]
#     print(len(conserved))
# NBS_A2B = [];NBS_B2A = [];NBS_conserved = []
# for i in NBS:
#     if i in A2B:
#         NBS_A2B.append(i)
#     elif i in B2A:
#         NBS_B2A.append(i)
#     elif i in conserved:
#         NBS_conserved.append(i)
# print(len(NBS_A2B),len(NBS_B2A),len(NBS_conserved),sep='\t')
# print('NBS A2B ratio is %.2f' % (len(NBS_A2B)/len(NBS)))
# print('NBS B2A ratio is %.2f' % (len(NBS_B2A)/len(NBS)))
# print('NBS Conserved ratio is %.2f' % (len(NBS_conserved)/len(NBS)))


#------------------2020.10.06 对表达量文件进行处理，求生物学重复的平均值------------------------------
# import pandas as pd
##-----------------AP85叶段数据-----------------------------------
# df = pd.read_excel(r'D:\Result\Transit\fpkm_caculate\全套\01.SES-208_2016.09-SingleEnd.fpkm_dataframe.xlsx')
# df ['Ss1'] = round((df['Ss1-1']+df['Ss2-1']+df['Ss3-1'])/3,2)
# df ['Ss2'] = round((df['Ss1-2']+df['Ss2-2']+df['Ss3-2'])/3,2)
# df ['Ss3'] = round((df['Ss1-3']+df['Ss2-3']+df['Ss3-3'])/3,2)
# df ['Ss4'] = round((df['Ss1-4']+df['Ss2-4']+df['Ss3-4'])/3,2)
# df ['Ss5'] = round((df['Ss1-5']+df['Ss2-5']+df['Ss3-5'])/3,2)
# df ['Ss6'] = round((df['Ss1-6']+df['Ss2-6']+df['Ss3-6'])/3,2)
# df ['Ss7'] = round((df['Ss1-7']+df['Ss2-7']+df['Ss3-7'])/3,2)
# df ['Ss8'] = round((df['Ss1-8']+df['Ss2-8']+df['Ss3-8'])/3,2)
# df ['Ss9'] = round((df['Ss1-9']+df['Ss2-9']+df['Ss3-9'])/3,2)
# df ['Ss10'] = round((df['Ss1-10']+df['Ss2-10']+df['Ss3-10'])/3,2)
# df ['Ss11'] = round((df['Ss1-11']+df['Ss2-11']+df['Ss3-11'])/3,2)
# df ['Ss12'] = round((df['Ss1-12']+df['Ss2-12']+df['Ss3-12'])/3,2)
# df ['Ss13'] = round((df['Ss1-13']+df['Ss2-13']+df['Ss3-13'])/3,2)
# df ['Ss14'] = round((df['Ss1-14']+df['Ss2-14']+df['Ss3-14'])/3,2)
# df ['Ss15'] = round((df['Ss1-15']+df['Ss2-15']+df['Ss3-15'])/3,2)
# df = df.iloc[:,[0,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60]]
# df.to_excel(r'D:\Result\Transit\fpkm_caculate\全套\SES_datnight_fpkm.xlsx',header=True,index=False)

##----------------AP85昼夜节律数据----------------------------------
# df = pd.read_excel(r'D:\Result\Transit\fpkm_caculate\全套\09.day_light_rhythm_2017.05.SES.fpkm_dataframe.xlsx')
###-------------2h--------------
# df['6:00'] =round((df['SE1-1_HJM2HALXX_L6']+df['SE1-2_HJML5ALXX_L8']+df['SE1-3_HJM2HALXX_L5'])/3,2)
# df['8:00'] =round((df['SE2-1_HGVFVALXX_L6']+df['SE2-2_HJM2HALXX_L5']+df['SE2-3_HJM2HALXX_L4'])/3,2)
# df['10:00'] =round((df['SE3-1_HJML5ALXX_L7']+df['SE3-2_HJM2HALXX_L4']+df['SE3-3_HJM2HALXX_L4'])/3,2)
# df['12:00'] =round((df['SE4-1_HJM2HALXX_L5']+df['SE4-2_HJML5ALXX_L7']+df['SE4-3_HJML5ALXX_L7'])/3,2)
# df['14:00'] =round((df['SE5-1_HH2MJALXX_L1']+df['SE5-2_HJM2HALXX_L1']+df['SE5-3_HJM2HALXX_L4'])/3,2)
# df['16:00'] =round((df['SE6-1_HJML5ALXX_L7']+df['SE6-2_HH2MJALXX_L1']+df['SE6-3_HJM2HALXX_L4'])/3,2)
# df['18:00'] =round((df['SE7-1_HH2MJALXX_L1']+df['SE7-2_HH2MJALXX_L1']+df['SE7-3_HJML5ALXX_L7'])/3,2)
# df['20:00'] =round((df['SE8-1_HH2MJALXX_L1']+df['SE8-2_HJM2HALXX_L6']+df['SE8-3_HH2MJALXX_L1'])/3,2)
# df['22:00'] =round((df['SE9-1_HJML5ALXX_L8']+df['SE9-2_HH2MJALXX_L1']+df['SE9-3_HJML5ALXX_L7'])/3,2)
# df['24:00'] =round((df['SE10-1_HJML5ALXX_L7']+df['SE10-2_HJM2HALXX_L7']+df['SE10-3_HH2MJALXX_L1'])/3,2)
# df['2:00s'] =round((df['SE11-1_HH2MJALXX_L1']+df['SE11-2_HJM2HALXX_L4']+df['SE11-3_HH2MJALXX_L1'])/3,2)
# df['4:00s'] =round((df['SE12-1_HH2MJALXX_L1']+df['SE12-2_HJML5ALXX_L8']+df['SE12-3_HJM2HALXX_L1'])/3,2)
# df = df.iloc[:,[0,60,61,62,63,64,65,66,67,68,69,70,71]]
# df.to_excel(r'D:\Result\Transit\fpkm_caculate\全套\SES_daynight_2h.xlsx',header=True,index=False)
###----------------------------
###--------------4h--------------
# df['6:00'] =round((df['S1-1_HGVJ2ALXX_L8']+df['S1-2_HGVJ2ALXX_L3']+df['S1-3_HGVJ2ALXX_L2'])/3,2)
# df['10:00'] =round((df['S2-1_HGVJ2ALXX_L3']+df['S2-2_HJLWYALXX_L2']+df['S2-3_HJLWYALXX_L2'])/3,2)
# df['14:00'] =round((df['S3-1_HGVJ2ALXX_L3']+df['S3-2']+df['S3-3_HGVJ2ALXX_L2'])/3,2)
# df['18:00'] =round((df['S4-1_HGVJ2ALXX_L1']+df['S4-2_HJM2NALXX_L1']+df['S4-3_HGVJ2ALXX_L7'])/3,2)
# df['22:00'] =round((df['S5-1_HJLWYALXX_L2']+df['S5-2_HGVJ2ALXX_L1']+df['S5-3_HJMTCALXX_L6'])/3,2)
# df['2:00s'] =round((df['S6-1_HGVJ2ALXX_L1']+df['S6-2_HGVJ2ALXX_L1']+df['S6-3_HGVJ2ALXX_L1'])/3,2)
# df['6:00s'] =round((df['S7-1_HGVJ2ALXX_L1']+df['S7-2_HGVJ2ALXX_L1']+df['S7-3_HJM2NALXX_L1'])/3,2)
# df = df.iloc[:,[0,60,61,62,63,64,65,66]]
# df.to_excel(r'D:\Result\Transit\fpkm_caculate\全套\SES_daynight_4h.xlsx',header=True,index=False)
# print(df)

##---------------NP-X-----------------------------------
# df = pd.read_excel(r'C:\Users\dell\Desktop\Nepal_project\Nepal_fpkm_df.raw.xlsx')
# df['L'] = round((df.iloc[:,1]+df.iloc[:,2]+df.iloc[:,3])/3,2)
# df['S'] = round((df.iloc[:,4]+df.iloc[:,5]+df.iloc[:,6])/3,2)
# df =df.iloc[:,[0,7,8]]
# df.to_excel(r'C:\Users\dell\Desktop\Nepal_project\Nepal_fpkm_df.xlsx',header=True,index=False)
# print(df)

##---------------AP85成熟期茎6和叶片-----------------------
# df = pd.read_excel(r'D:\Result\Transit\fpkm_caculate\全套\01.SES-208_2016.09-PairedEnd.fpkm_dataframe.xlsx')
# df['L'] = df['m-SES-mature-leaf-3']
# df['S6'] = df['m-SES-stem6-1']
# df = df.iloc[:,[0,50,51]]
# df.to_excel(r'C:\Users\dell\Desktop\Nepal_project\20201007_AP85_NPX_fpkm\AP85_fpkm.xlsx',header=True,index=False)
# print(df)

#---------------制作AP852019版本的GOannotation------------

# Tair_df = pd.read_table(r'C:\Users\dell\Desktop\脚本测试\20200819 clusterProfiler + eggNOG\Sspon_all_genome\ATH_GO_GOSLIM.txt',sep='\t',skiprows=4,header=None)
# df = Tair_df[Tair_df.iloc[:,1].str.contains('gene')]
# df = df.iloc[:,[2,5]]
# df.rename(columns={2:'gene_id',5:'GO'},inplace=True)
# # print(df)
# Dict = df.set_index('gene_id').to_dict()['GO']
# print(df)
# AP85_anno = pd.read_table(r'C:\Users\dell\Desktop\脚本测试\20200819 clusterProfiler + eggNOG\Sspon_all_genome\SsponAnno_replaced.v20190122.txt',sep='\t')
# AP85_anno['GOid'] = AP85_anno['Best-hir-arabi-name'].map(Dict)
# AP85_anno_df = AP85_anno.iloc[:,[0,-1,4]]
# print(AP85_anno_df)
# AP85_anno_df = AP85_anno_df.dropna(how='any')
# AP85_anno_df.to_csv(r'C:\Users\dell\Desktop\脚本测试\20200819 clusterProfiler + eggNOG\Sspon_all_genome\AP85_19version_GOanno.tab',sep='\t',index=False)
# # print(df)


#--------------concat--------------------------------
import os
# os.chdir(r'C:\Users\dell\Desktop\Nepal_project\20201007_AP85_NPX_fpkm')
# file_list = os.listdir('./')
# print(file_list[1:-2])
# # file1 = pd.read_table(file_list[2],comment='#')
# # print(file1)
# concat_list = []
#
# for file in file_list[1:-2]:
#     df = pd.read_table(file,sep='\t',comment='#',header=None)
#     print(df)
#     concat_list.append(df)
#
#
# result = pd.concat(concat_list)
# print(result)
# result.to_csv(r'C:\Users\dell\Desktop\Nepal_project\20201007_AP85_NPX_fpkm\Nepal_v_AP85_concat.anchors',sep='\t',index=False,header=None)


#--------------2020.10.13 renamePhotosynthesisGene----------------------------

# rename_df = pd.read_table(r'D:\调控元件生信任务\数据\rename_info.list',names=['v18','v19'])
# print(rename_df)
# Dict = rename_df.set_index('v18').to_dict()['v19']
# print(Dict)
# PhotosynthesisGene = pd.read_table(r'D:\调控元件生信任务\数据\Photosynthesis_gene_18v.txt',names=['v18'])
# print(PhotosynthesisGene)
# PhotosynthesisGene['v19'] = PhotosynthesisGene['v18'].map(Dict)
# print(PhotosynthesisGene)
# # PhotosynthesisGene.iloc[:,1].to_csv(r'D:\调控元件生信任务\数据\Photosynthesis_gene_19v.txt',header=None,index=False)
# query_gene_list = PhotosynthesisGene['v19'].tolist()
# with open(r'D:\调控元件生信任务\数据\Photosynthesis_gene_19v.txt','w') as f1:
#     for i in query_gene_list:
#         if str(i) != 'nan' :
#             print(i)
#             f1.write(i+'\n')


###

import os
# pd.set_option('display.max_columns',12)
#----------------------筛选未鉴定到调控元件的基因ID-----------------------
# Sspon_all_gene = pd.read_table(r'C:\Users\admin\Desktop\genome-wide_geneID_v0806.txt',header=None)
# list1 = Sspon_all_gene.iloc[:,0].tolist()
# #
# screende = []
# with open(r'C:\Users\admin\Desktop\screened_Sspon_geneV2.txt') as f1:
#     for line in f1:
#         line = line.strip()
#         if not line.startswith('GeneID'):
#             if line.endswith('+'):
#                 screende.append(line[:-1])
#             elif line.endswith('-'):
#                 screende.append(line[:-1])
#             elif not line.endswith('+'):
#                 if not line.endswith('-'):
#                     screende.append(line)
#
# for gene in screende:
#     if not gene in list1:
#         print(gene)
        #
        #     if not line.endswith('+'):
        #         if not line.endswith('-'):
        #             print(line)
# -------------------------------------------------------------------

#---------------------将未鉴定元件的基因的plantcare输出结果进行清洗-------------
# df = pd.read_table(r'D:\chrom download\PlantCARE_12559_unscreen_plantCARE\plantCARE_output_PlantCARE_12559.tab',sep='\t',index_col=None,header=None)
# # df = df.reset_index()
# df.rename(columns={0: 'GeneID', 1: 'Site Name', 2: 'Sequence', 3: 'Position', 4: 'Length', 5: 'Strand', 6: 'Organism',
#                  7: 'Function'}, inplace=True)
# df = df.dropna(how='any')
# df = df.groupby(['GeneID', 'Site Name', 'Sequence', 'Strand', 'Length', 'Organism', 'Function'])['Position'].apply(
#         list)
# df = df.reset_index()
#
# print(df)
#-----------------------------------------------------------------

#---------------------------Bru1 blast结果分析----------------------------------
# df = pd.read_table(r'C:\Users\admin\Desktop\硕士课题\Bru1\Bru1_LA.out',comment='#',names=['query_acc','subject_acc','identity','alignment_length',
#                                                                                    'mismatches','gap_opens','q_start','q_end','s_start','s_end',
#                                                                                    'evalue','bit_score'],float_precision='round_trip')
# print(df)
# df = df.sort_values(['bit_score'],ascending=False)
# df = df[df['evalue'] == 0]
# print(df)
# df = df.to_csv(r'C:\Users\admin\Desktop\硕士课题\Bru1\Bru1_LA_clean.out',sep='\t',header=True,index=False)
#-----------------------------------------------------------------------------

# --------------------------------20201105给老张做ppt--------------------------------------------------
# import re
# df = pd.read_table(r'E:\潘浩然\调控元件生信任务\数据\Sspon_all_genome_CREv1101.tab',sep='\t')
# # print(df)
# df['GeneID'] = df['GeneID'].str[:-1]
# # print(df)
# photosynthesis = pd.read_table(r'D:\pycharm\panhr\FAFU_CGB\Genomic_exercise\CRE_biotask\allele4_gene.txt',sep='\t',header=None)
# # print(photosynthesis.iloc[:,0].tolist())
# # print(len(photosynthesis.iloc[:,0].tolist()))
#
#
# df1 = df[df['GeneID'].isin(photosynthesis.iloc[:,0].tolist())]
# print(df1)
# # df1.to_csv(r'D:\pycharm\panhr\FAFU_CGB\Genomic_exercise\CRE_biotask\allele3_gene.tab',sep='\t',header=True,index=False)
# print(len(df1.iloc[:,0].unique()))
# #
# df1 = df1.iloc[:,[0,1,-1]]
# df1['Count'] = df1['Position'].str.count(',')+1
# print(df1)
# df1 = df1.groupby(['GeneID','Site Name'])['Count'].apply(list)
# df1 = df1.reset_index()
# print(df1)   #TODO 将count中的数值相加在新的一列中显示
#
# # df1.to_excel(r'D:\pycharm\panhr\FAFU_CGB\Genomic_exercise\CRE_biotask\allele_gene4.xlsx',header=True,index=False)
# 使用EXCEL进行处理后在读取
# df2 = pd.read_excel(r'D:\pycharm\panhr\FAFU_CGB\Genomic_exercise\CRE_biotask\allele_gene4.xlsx')
# print(df2)
# df2 = df2.groupby(by='Site Name').agg({"Count":sum}).reset_index()
# df2 = df2.sort_values(by='Count',ascending=False)
# print(df2)
#
# df2.to_csv(r'D:\pycharm\panhr\FAFU_CGB\Genomic_exercise\CRE_biotask\allele4_result.tab',sep='\t',header=True,index=False)
#---------------------------------------------------------------------------------------------

##----------------LA昼夜节律数据----------------------------------
# df = pd.read_excel(r'E:\潘浩然\Result\Transit\fpkm_caculate\全套\09.day_light_rhythm_2017.05.LA.fpkm_dataframe.xlsx')
##-------------2h--------------
# df['6:00'] =round((df['SE1-1_HJM2HALXX_L6']+df['SE1-2_HJML5ALXX_L8']+df['SE1-3_HJM2HALXX_L5'])/3,2)
# df['8:00'] =round((df['SE2-1_HGVFVALXX_L6']+df['SE2-2_HJM2HALXX_L5']+df['SE2-3_HJM2HALXX_L4'])/3,2)
# df['10:00'] =round((df['SE3-1_HJML5ALXX_L7']+df['SE3-2_HJM2HALXX_L4']+df['SE3-3_HJM2HALXX_L4'])/3,2)
# df['12:00'] =round((df['SE4-1_HJM2HALXX_L5']+df['SE4-2_HJML5ALXX_L7']+df['SE4-3_HJML5ALXX_L7'])/3,2)
# df['14:00'] =round((df['SE5-1_HH2MJALXX_L1']+df['SE5-2_HJM2HALXX_L1']+df['SE5-3_HJM2HALXX_L4'])/3,2)
# df['16:00'] =round((df['SE6-1_HJML5ALXX_L7']+df['SE6-2_HH2MJALXX_L1']+df['SE6-3_HJM2HALXX_L4'])/3,2)
# df['18:00'] =round((df['SE7-1_HH2MJALXX_L1']+df['SE7-2_HH2MJALXX_L1']+df['SE7-3_HJML5ALXX_L7'])/3,2)
# df['20:00'] =round((df['SE8-1_HH2MJALXX_L1']+df['SE8-2_HJM2HALXX_L6']+df['SE8-3_HH2MJALXX_L1'])/3,2)
# df['22:00'] =round((df['SE9-1_HJML5ALXX_L8']+df['SE9-2_HH2MJALXX_L1']+df['SE9-3_HJML5ALXX_L7'])/3,2)
# df['24:00'] =round((df['SE10-1_HJML5ALXX_L7']+df['SE10-2_HJM2HALXX_L7']+df['SE10-3_HH2MJALXX_L1'])/3,2)
# df['2:00s'] =round((df['SE11-1_HH2MJALXX_L1']+df['SE11-2_HJM2HALXX_L4']+df['SE11-3_HH2MJALXX_L1'])/3,2)
# df['4:00s'] =round((df['SE12-1_HH2MJALXX_L1']+df['SE12-2_HJML5ALXX_L8']+df['SE12-3_HJM2HALXX_L1'])/3,2)
# df = df.iloc[:,[0,60,61,62,63,64,65,66,67,68,69,70,71]]
# df.to_excel(r'D:\Result\Transit\fpkm_caculate\全套\SES_daynight_2h.xlsx',header=True,index=False)
# print(df)
# print(df)
# df['6:00'] = df.iloc
###----------------------------
###--------------4h--------------
# df['6:00'] =round((df['S1-1_HGVJ2ALXX_L8']+df['S1-2_HGVJ2ALXX_L3']+df['S1-3_HGVJ2ALXX_L2'])/3,2)
# df['10:00'] =round((df['S2-1_HGVJ2ALXX_L3']+df['S2-2_HJLWYALXX_L2']+df['S2-3_HJLWYALXX_L2'])/3,2)
# df['14:00'] =round((df['S3-1_HGVJ2ALXX_L3']+df['S3-2']+df['S3-3_HGVJ2ALXX_L2'])/3,2)
# df['18:00'] =round((df['S4-1_HGVJ2ALXX_L1']+df['S4-2_HJM2NALXX_L1']+df['S4-3_HGVJ2ALXX_L7'])/3,2)
# df['22:00'] =round((df['S5-1_HJLWYALXX_L2']+df['S5-2_HGVJ2ALXX_L1']+df['S5-3_HJMTCALXX_L6'])/3,2)
# df['2:00s'] =round((df['S6-1_HGVJ2ALXX_L1']+df['S6-2_HGVJ2ALXX_L1']+df['S6-3_HGVJ2ALXX_L1'])/3,2)
# df['6:00s'] =round((df['S7-1_HGVJ2ALXX_L1']+df['S7-2_HGVJ2ALXX_L1']+df['S7-3_HJM2NALXX_L1'])/3,2)
# df = df.iloc[:,[0,60,61,62,63,64,65,66]]
# df.to_excel(r'D:\Result\Transit\fpkm_caculate\全套\SES_daynight_4h.xlsx',header=True,index=False)
# print(df)

# def silding_windows_bed(bedfile):
#     Chrdict = {}
#     site_num = 0
#     with open(bedfile) as f:
#         for line in f:
#             line_list = line.strip().split('\t')
#             if line_list[0] not in Chrdict:
#                 Chrdict[line_list[0]] = 0
#                 site_num +=1
#             elif line_list[0] in Chrdict:
#
# import sys
# def count_CG_in_fasta(fasta_file):
#     GC_num = 0
#     with open(fasta_file) as f:
#         for line in f:
#             if line.startswith('>'):
#                 continue
#             else:
#                 GC_num += line.count('GC')
#                 GC_num += line.count('CG')
#     print(GC_num)
# count_CG_in_fasta(sys.argv[1])
from path import Path
import sys

def find_SCP_geneID(file_path):
    file_list = Path(file_path).files()
    for file in file_list:
        with open(file) as f:
            lines = (line.strip() for line in f)
            for line in lines:
                if line.startswith('>'):
                    print(line.strip()[1:])
        # print('-'*20)

# if __name__ == '__main__':
#     if len(sys.argv) !=2:
#         print('Usage:')
#         print('\tpython {0} SCPfa_path > output.txt'.format(sys.argv[0]))
#     else:
#         find_SCP_geneID(sys.argv[1])
import argparse
def Extract_promoter_Seq(species_ID,promoter_file,output):
    fa_dict = {}
    strand_dict = {}
    with open(promoter_file) as f:
        lines = (line.strip() for line in f)
        for line in lines:
            if line.startswith('>'):
                key = line[1:].split()[0]
                fa_dict[key] = ''
                strand_dict[key] = line[1:].split()[1]
            else:
                fa_dict[key] = line
    with open(output,'w') as w:
        for gene in species_ID:
            if gene in fa_dict.keys():
                w.write('>'+gene+strand_dict[gene]+'\n')
                w.write(fa_dict[gene]+'\n')
            else:
                continue


def clean_geneID(geneID_file):
    specie1_geneID = [];specie2_geneID = [];specie3_geneID = []
    with open(geneID_file) as f:
        lines = (line.strip() for line in f)
        for line in lines:
            if line.startswith('LOC_Os'):
                line = line.split('.')[0]
                specie1_geneID.append(line)
            elif line.startswith('Sobic'):
                line = '.'.join(line.split('.')[0:2])
                specie2_geneID.append(line)
            elif line.startswith('Zm'):
                line = line.split('_')[0]
                specie3_geneID.append(line)
    return specie1_geneID, specie2_geneID, specie3_geneID

# if __name__ == '__main__':
#     if len(sys.argv) != 5:
#         print('Find the promoter of the query gene(3 species)')
#         print('Usage:')
#         print('\tpython {0} geneID_file 1/2/3 promoter_file output_file')
#     else:
#         species1ID, species2ID, species3ID = clean_geneID(sys.argv[1])
#         if sys.argv[2] == '1':
#             Extract_promoter_Seq(species1ID,sys.argv[3],sys.argv[4])
#         if sys.argv[2] == '2':
#             Extract_promoter_Seq(species2ID,sys.argv[3],sys.argv[4])
#         if sys.argv[2] == '3':
#             Extract_promoter_Seq(species3ID,sys.argv[3],sys.argv[4])

def find_promoter_pos(geneposbed):
    '''
    geneposbed 5 col : Chr,start,end,strand,geneid
    :param geneposbed:
    :return:
    '''
    with open(geneposbed) as f:
        for line in f :
            line_list = line.strip().split()
            if line_list[3] == '+':
                print(line_list[0],int(line_list[1])-2000,line_list[1],line_list[3],line_list[4],sep='\t')
            if line_list[3] == '-':
                print(line_list[0],line_list[2],int(line_list[2])+2000,line_list[3],line_list[4],sep='\t')

# if __name__ == '__main__':
#     if len(sys.argv) != 2:
#         print('Usage:')
#         print('\tpython {0} GenePosBed > PromoterPosBed'.format(sys.argv[0]))
#     else:
#         find_promoter_pos(sys.argv[1])
# find_promoter_pos(r'C:\Users\admin\Desktop\Chr3A.gene.bed')
def clean_CpGGene_interserct(inter_file):
    line_dict = {}
    with open(inter_file) as f:
        for line in f:
            line_list = line.strip().split('\t')
            key = line_list[4]
            if not key in line_dict:
                line_dict[key] = int(line_list[8])
            else:
                line_dict[key] += int(line_list[8])

            # print(key)
    for i in line_dict.keys():
        print(i,line_dict[i],sep='\t')


# clean_CpGGene_interserct(r'C:\Users\admin\Desktop\20201205\CGI\allgenome.inter.bed')

from rich import print
from rich.console import Console
from rich import pretty

print("Hello, [bold magenta]World[/bold magenta]!", ":vampire:")
console = Console()
console.print('hello')
console.print("Hello", "World!", style="bold red")
console.print("Hello World!", style="bold red")

from rich.console import Console
from rich.table import Column, Table

console = Console()

table = Table(show_header=True, header_style="bold magenta")
table.add_column("Date", style="dim", width=12)
table.add_column("Title")
table.add_column("Production Budget", justify="right")
table.add_column("Box Office", justify="right")
table.add_row(
    "Dev 20, 2019", "Star Wars: The Rise of Skywalker", "$275,000,000", "$375,126,118"
)
table.add_row(
    "May 25, 2018",
    "[red]Solo[/red]: A Star Wars Story",
    "$275,000,000",
    "$393,151,347",
)
table.add_row(
    "Dec 15, 2017",
    "Star Wars Ep. VIII: The Last Jedi",
    "$262,000,000",
    "[bold]$1,332,539,889[/bold]",
)

console.print(table)

from rich.progress import track
import time
# for step in track(range(100)):
#     # print(step)
#     time.sleep(0.2)
from rich.console import Console
from rich.syntax import Syntax

my_code = '''
def iter_first_last(values: Iterable[T]) -> Iterable[Tuple[bool, bool, T]]:
    """Iterate and generate a tuple with a flag for first and last value."""
    iter_values = iter(values)
    try:
        previous_value = next(iter_values)
    except StopIteration:
        return
    first = True
    for value in iter_values:
        yield first, False, previous_value
        first = False
        previous_value = value
    yield first, True, previous_value
'''
syntax = Syntax(my_code, "python", theme="monokai", line_numbers=True)
console = Console()
console.print(syntax)

