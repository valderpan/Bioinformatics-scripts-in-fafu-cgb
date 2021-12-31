#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/1/23


from path import Path
from collections import Counter
import pandas as pd
import re



def get_chr_base_matirx(matirx_file,start,end,output):
    df = pd.read_table(matirx_file,header=None,sep='\t',names=['firstbin','secondbin','frequency'],float_precision='round_trip')
    df_start = df[df['firstbin'] == int(start)].index[0]
    df_end = df[df['firstbin'] == int(end)].index[-1]
    df = df.iloc[df_start:df_end+1,]
    print(df)
    df.to_csv(output,sep='\t',header=True,index=False)


# get_chr_base_matirx(r'E:\Nepal_project\20210123_Chr_recombination\Nepal_500000_iced.matrix')
get_chr_base_matirx(r'E:\Nepal_project\20210123_Chr_recombination\Nepal_500000_iced.matrix','3600','3673',
                    r'E:\Nepal_project\20210123_Chr_recombination\Nepal_Chr5D_Recomd_6AHiC.txt')
from path import Path
import pandas as pd
files = Path(r'E:\Badila\synteny_analysis\C03').files()
queryfiles = [file for file in files if 'nonconserverd_group' in file]
df1 = pd.read_excel(r'E:\Badila\synteny_analysis\C03\C03.sort.xlsx',header=None)
print(df1.shape[0])
# print(queryfiles)
num = 0
for file in queryfiles:

    for i in range(0,20):
        try:
            df = pd.read_excel(file,sheet_name=i,header=None)
            # print(df)
            num += df.shape[0]
        except:
            continue
print(num)

for i in range(0,20):
    try:
        df2 = pd.read_excel(r'E:\Badila\synteny_analysis\C03\conserverd_group.xlsx',sheet_name=i,header=None)
        # print(df2)
        num += df2.shape[0]
    except:
        continue
print(num)



import pandas as pd

df1 = pd.read_excel(r'E:\Badila\synteny_analysis\C03\conserverd_group.xlsx',header=None,sheet_name=3)
print(df1)
df1_gene = df1.iloc[:,0].tolist()
df2 = pd.read_excel(r'E:\Badila\synteny_analysis\C03\C03.sort.tmp.xlsx',header=None,sheet_name=1)
print(df2)
df2_gene = df2.iloc[:,0].tolist()

for i in df2_gene:
    if not i in df1_gene:
        print(i)

import pandas as pd
#
df = pd.read_table(r'E:\潘浩然\调控元件生信任务\analysis\callpeaks\75150\SES208_peaks.narrowPeak.txt',sep='\t',header=None,float_precision='round_trip')
# df = df[~df.iloc[:,0].str.contains('tig')]
# print(df)
df.rename(columns={0:'ChrID',1:'idr_start',2:'idr_end'},inplace=True)
# print(df)
df = df.sort_values(by=['ChrID','idr_start'],ignore_index=True)
# print(df)

df =df.iloc[:,[0,1,2]]
peak_name = ['SES208_peak_{}'.format(i) for i in range(1,df.shape[0]+1)]

df['peak_name'] = peak_name

print(df)
df.to_csv(r'E:\潘浩然\调控元件生信任务\analysis\callpeaks\75150\SES208_peaks.narrowPeak.idr.txt',sep='\t',header=False,index=False)

#--------统计高粱基因组每条染色体上有多少基因-------------


#--------------tab2xlsx-------------------------
from path import Path
import re
files = [file for file in Path('E:\Badila\synteny_analysis\解麻子').files() if file.endswith('.info')]
# print(files)
for file in files:
    df = pd.read_table(file,sep='\t',names=['V1','V2','V3','V4','V5','V6','V7'])
    file_name = re.findall('([0-9a-zA-Z\.\-]+)\.info',file)[0]
    # print(df)
    # print(file_name)
    # break
    df.to_excel(file_name+'.xlsx',header=False,index=False)
print('start')

import sys
Dict = {}
with open(sys.argv[1]) as f:
    lines = (line.strip() for line in f)
    for line in lines:
        if line.startswith('S'):
            Dict[line.split('\t')[1]] = line.split('\t')[2]
        else:
            Dict[line.split('\t')[1]] += line
for key in Dict.keys():
    print('>'+key)
    print(Dict[key])





a='ChrSy.fgenesh.mRNA.74'
print(re.sub('mRNA','gene',a))
print(len([dir for dir in Path('./').dirs() if dir.basename() == 'cds_SCG']))
print(Path('./').dirs())
df = pd.read_table(r'E:\潘浩然\调控元件生信任务\analysis\Sugarcane\OsZmSb\Maize_7days_leaf_ACRs.bed',header=None)
print(df)
df_chr = df[~df[0].str.contains('B73')]
print(df_chr)
df_chr[0] = 'Chr'+df_chr[0]
print(df_chr)
result = pd.concat([df_chr,df[df[0].str.contains('B73')]])
print(result)
result.to_csv(r'E:\潘浩然\调控元件生信任务\analysis\Sugarcane\OsZmSb\Maize_7days_leaf_ACRs.corrected.bed',index=False,header=False,sep='\t')


from Bio import SeqIO
import sys
def Fa2dict_Bio(fasta_file):
    fasta_dict = {}
    for seq in SeqIO.parse(fasta_file,'fasta'):
        fasta_dict[seq.id] = seq.seq
    return fasta_dict

def get_query_ID(ID_file):
    ID_list = []
    with open(ID_file,'r') as f:
        lines = (line.strip() for line in f)
        for line in lines:
           ID_list.append(line)
    print('Query_ID numbers : {}'.format(len(ID_list)))
    return ID_list

def get_query_seq(ID_list,fasta_Dict):
    for id in ID_list:
        if id in fasta_Dict.keys():
            print('>'+id)
            print(str(fasta_Dict[id]))

if __name__ == '__main__':
    D = Fa2dict_Bio(sys.argv[1])
    L = get_query_ID(sys.argv[2])
    get_query_seq(L,D)


def parse_kemerfreq(file):
    kmer_dict = {}
    total_lines = 0
    with open(file,'r') as f:
        lines = (line.strip() for line in f)
        for line in lines:
            # kmer.append(line)
            kmer,freq = line.split('\t')
            kmer_dict[int(kmer)] = int(freq)
            total_lines+=1
    # print(total_lines)
    for i,j in kmer_dict.items():
        # print(type(i),type(j))
        if int(i)+1 < total_lines and i -1 >0:
            # if int(kmer_dict[str(int(i)+1)]) < int(kmer_dict[i]) and int(kmer_dict[str(int(i)-1)]) < int(kmer_dict[i]):
            if kmer_dict[i+1]< kmer_dict[i] and kmer_dict[i-1]< kmer_dict[i]:
                if kmer_dict[i+2]< kmer_dict[i] and kmer_dict[i-2]< kmer_dict[i]:
                    if kmer_dict[i + 3] < kmer_dict[i] and kmer_dict[i - 3] < kmer_dict[i]:
                        if kmer_dict[i + 4] < kmer_dict[i] and kmer_dict[i - 4] < kmer_dict[i]:
                            if kmer_dict[i + 5] < kmer_dict[i] and kmer_dict[i - 5] < kmer_dict[i]:
                                if kmer_dict[i + 6] < kmer_dict[i] and kmer_dict[i - 6] < kmer_dict[i]:
                                    if kmer_dict[i + 7] < kmer_dict[i] and kmer_dict[i - 7] < kmer_dict[i]:
                                        if kmer_dict[i + 8] < kmer_dict[i] and kmer_dict[i - 8] < kmer_dict[i]:
                                            if kmer_dict[i + 9] < kmer_dict[i] and kmer_dict[i - 9] < kmer_dict[i]:
                                                if kmer_dict[i + 10] < kmer_dict[i] and kmer_dict[i - 10] < kmer_dict[i]:
                                                    if i <= 1000:
                                                        print(i)
            # break
parse_kemerfreq(r'E:\中转站\freq.stat.2colum')
from collections import Counter


def Chr2stablegroup(refChrNum):
    Chr2groupDict = {}
    for chr in range(1,refChrNum+1):
        if chr == 1:
            Chr2groupDict['Chr{}'.format(chr)] = ['group'+ str(x) for x in [4,7,14,31,37,47,52,65]]
        elif chr == 2:
            Chr2groupDict['Chr{}'.format(chr)] = ['group'+ str(x) for x in [3,6,15,20,21,25,54,63]]
        elif chr == 3:
            Chr2groupDict['Chr{}'.format(chr)] = ['group'+ str(x) for x in [1,11,12,19,22,24,59,73]]
        elif chr == 4:
            Chr2groupDict['Chr{}'.format(chr)] = ['group'+ str(x) for x in [17,23,26,35,40,57,62,69]]
        elif chr == 5:
            Chr2groupDict['Chr{}'.format(chr)] = ['group'+ str(x) for x in [36,39,43,46,50,68,71,77]]
        elif chr == 6:
            Chr2groupDict['Chr{}'.format(chr)] = ['group'+ str(x) for x in [5,18,44,51,61,67,70,72]]
        elif chr == 7:
            Chr2groupDict['Chr{}'.format(chr)] = ['group'+ str(x) for x in [9,10,28,33,34,66,79,80]]
        elif chr == 8:
            Chr2groupDict['Chr{}'.format(chr)] = ['group'+ str(x) for x in [2,27,49,53,55,56,75,78]]
        elif chr == 9:
            Chr2groupDict['Chr{}'.format(chr)] = ['group'+ str(x) for x in [8,13,30,32,38,60,64,74]]
        elif chr == 10:
            Chr2groupDict['Chr{}'.format(chr)] = ['group'+ str(x) for x in [16,29,41,42,45,48,58,76]]
    return Chr2groupDict

def stat_assembly_chrGeneNum(group_gff3,D):
    df = pd.read_table(group_gff3,comment='#',header=None)
    gene_df = df[df[2] == 'gene']
    print(gene_df)
    gene_df_num = gene_df[0].tolist()
    print(len(gene_df_num))
    b = dict(Counter(gene_df_num))
    print(b)
    bb = pd.DataFrame([b]).T
    bb = bb.reset_index()
    print(bb)
    type_list = []
    for index,row in bb.iterrows():
        # print(row['index'])
        for chr in D.keys():
            if row['index'] in D[chr]:
                type_list.append(chr)
    print(len(type_list))
    bb['type'] = type_list
    print(bb)
    bb.to_excel(r'E:\Badila\synteny_analysis\Badila染色体第一次校正\矫正后\each_Sbchr_geneNumber.xlsx',header=False,index=False)

if __name__ == '__main__':
    D = Chr2stablegroup(10)
    print(D)
    stat_assembly_chrGeneNum(r'E:\中转站\Badila.group.gff3',D)

df_origin = pd.read_excel(r'D:\pycharm\panhr\FAFU_CGB\Genomic_exercise\ROC22\CPCS\group73.origin_df.xlsx',header=None)
df_remove = pd.read_excel(r'D:\pycharm\panhr\FAFU_CGB\Genomic_exercise\ROC22\CPCS\group73.remove_df.xlsx',header=None)
df_removed = pd.read_excel(r'D:\pycharm\panhr\FAFU_CGB\Genomic_exercise\ROC22\CPCS\group73.removed_df.xlsx',header=None)
df_add = pd.read_excel(r'D:\pycharm\panhr\FAFU_CGB\Genomic_exercise\ROC22\CPCS\group73.add_df.xlsx',header=None)
print(df_remove)
print(df_removed)
print(df_add)
print(df_origin)
# for i in df_remove[0].tolist():
#     if i in df_removed[0].tolist():
#         print(i)
df_remove_ctg = df_remove[0].tolist()
print(len(df_remove_ctg))

a = df_origin[df_origin[0].isin(df_remove_ctg)]
print(a.shape[0])

for i in df_remove_ctg:
    if i not in df_origin[0].tolist():
        print(i)


from BCBio import GFF
from collections import Counter
biotype = []
in_handle = open(r'E:\潘浩然\调控元件生信任务\analysis\callpeak\Sspon.v20190103.gff3','r')
geneID2loc = {}
num1 = 0;num2 = 0
for rec in GFF.parse(in_handle):
    # print(rec)
    num1 +=1
    for feature in rec.features:
        # print(feature)
        # print(feature.qualifiers)
        num2+=1
        for subfeature in feature.sub_features:
            print(subfeature)
        # print(feature.location)
        # print(type(feature.location))
        # geneID2loc[feature.id] = str(feature.location)
    # break
        # if feature.type == "gene":
        #     # print(feature.qualifiers) # {'ID': ['Sspon.08G0030880-1D'], 'Name': ['Sspon.08G0030880-1D'], 'source': ['FGENESH']}
        #     biotype.append(feature.qualifiers['source'][0])
        break
in_handle.close()
print(num1)
print(num2)
# print(len(biotype))
# print(biotype[0:15])
# print(Counter(biotype))

def merge_ks_result(file_path):
    files = [file for file in Path(file_path).files() if file.endswith('.result')]
    writer = pd.ExcelWriter('LA_allele_Ks.xlsx')
    for file in files:
        file_prefix = re.findall('\_([A-Z\_]+)\.result',file.basename())[0]
        df = pd.read_table(file,sep='\t')
        df = df.iloc[:,[0,3]]
        df.to_excel(writer,sheet_name=file_prefix,index=False)
        # print(file_prefix)
    writer.save()
merge_ks_result(r'E:\Nepal_project\20210616_NpX_review\LA_allele_Ks')

import sys
from path import Path

def stat_dup_ctgIngroup(file_path):
    files = [file for file in Path(file_path).files() if file_path.endswith('.txt')]
    concat_list = []
    for file in files:
        df = pd.read_table(file)
        # df = df.
        concat_list.append(df)
    concat_df = pd.concat(concat_list)
    return concat_df

def check_dup(df):
    df_dup = df['#Contig'].duplicated()
    print(df_dup)

if __name__ == '__main__':
    stat_dup_ctgIngroup(sys.argv[1])

def find_allele4_60DEG_ACRs(peak_file,DEG60geneID):
    geneID = []
    with open(DEG60geneID,'r') as f:
        lines = (line.strip() for line in f)
        for line in lines:
            geneID.append(line)
    print(geneID)

    df = pd.read_excel(peak_file)
    gene_prefix = df['geneId'].str.split('-',expand=True).iloc[:,0]
    df['gene_prefix'] = gene_prefix
    print(df)

    query_df = df[df['gene_prefix'].isin(geneID)]
    print(query_df)
    print(len(query_df['geneId'].unique().tolist()))
    a = query_df.groupby('category').agg('count')
    print(a)
    # query_df.to_excel('DEG60.peak.anno.xlsx',header=True,index=False)

    DEG60_noACRs = []
    for i in geneID:
        if i not in df['gene_prefix'].tolist():
            DEG60_noACRs.append(i)
    print(DEG60_noACRs)

if __name__ == '__main__':
    find_allele4_60DEG_ACRs(r'E:\潘浩然\调控元件生信任务\analysis\Sugarcane\peaks\SES208_peak_category.xlsx',r'E:\潘浩然\调控元件生信任务\analysis\Sugarcane\peaks\allele\DEG\allele4jvenn_geneID.txt')

import sys
def extract_gene_pos(gff_file):
    with open(gff_file,'r') as f:
        lines = (line.strip() for line in f )
        for line in f:
            if not line.startswith('#'):
                line_list = line.split()
                if line_list[2] == 'mRNA':
                    gene_ID = line_list[8].split(';')[1].split('=')[1]
                    print(line_list[0],line_list[2],line_list[3],line_list[4],gene_ID,sep='\t',)
if __name__ == '__main__':
    extract_gene_pos(sys.argv[1])
import math
import numpy as np
df = pd.read_table(r'C:\Users\admin\Desktop\NpX.v0804.chrscale_meanFPKM_Leaf.bed',sep='\t',names=['seqid','start','end','L','S'])
df['start'] = df['start']+1
df_L = df.iloc[:,[0,1,2,3]]
df_L['LogL'] = round(df_L['L'].apply(np.log10),3)
df_S = df.iloc[:,[0,1,2,4]]
df_S['LogS'] = round(df_S['S'].apply(np.log10),3)
print(df)
print(df_L)
print(df_S)
df_L.to_csv('NpX.v0804.chrscale_meanFPKM_Leaf.bed',sep='\t',header=False,index=False)
df_S.to_csv('NpX.v0804.chrscale_meanFPKM_Stem.bed',sep='\t',header=False,index=False)



import math
import numpy as np
import sys
import pandas as pd
df = pd.read_table(sys.argv[1],names=['seqid','start','end','num'])
df['lognum'] = round(df['num'].apply(np.log10),3)
print(df)
df.loc[:,['seqid','start','end','lognum']].to_csv(sys.argv[2],sep='\t',header=False,index=False)


df = pd.read_table(sys.argv[1],comment='#',names=['seqid','start','end','count'])
df['new'] = df['count']*100
df.loc[:,['seqid','start','end','new']].to_csv(sys.argv[2],sep='\t',header=False,index=False)

from Bio.Seq import Seq
a = Seq('abcdesf')
print(str(a))
print(str(a)[0:4])

df = pd.read_excel(r'E:\潘浩然\调控元件生信任务\analysis\Sugarcane\Ss_LA\LA-pruple_peakAnno.mutate.3ACR.xlsx')
uniq_gene = df['geneId'].unique().tolist()
print(len(uniq_gene))
len = df['width'].agg([sum])
print(len)

import sys
import pandas as pd
df = pd.read_excel(sys.argv[1])
# for i in range(1,10):
#     for j in ['A','B','C','D']:
#         qdf = df[df['seqnames'] == 'Chr{}{}'.format(i,j)]
#         qdf_a = qdf[qdf['score'] >=0]
#         qdf_b = qdf[qdf['score'] < 0]
df['width'] = df['end'] - df['start']
df = df.loc[:,['seqnames','start','end','width','strand','score']]
df.to_csv('NpX_v0630.allchr.flipped.compartment.bed',sep='\t',header=True,index=False)

def cal_density(file):
    df = pd.read_table(file)
    df = df[df['MyA'] > 0]
    for i in df['species'].unique().tolist():
        qdf = df[df['species'] == i]
        qdf_d = Counter(qdf['MyA'])
        print(i,qdf_d.most_common(10))
        if i == 'Af':
            print(qdf_d)

    # print(df)

cal_density(r'E:\FLL\Repeat_Sequence\LTR\FLL.LTR_insertTime.bed')


import sys
from path import Path
import pandas as pd
import re
files = [i for i in Path(sys.argv[1]).files() if i.endswith('.tab')]
I = []
for file in files:
    n = re.findall('([0-9a-zA-Z\-]+)\.tab',file)[0]
    df = pd.read_table(file,sep='\t')
    df[n] = df['FPKM']
    I.append(df[n])

res = pd.concat(I,axis=1)
print(res)
res.to_csv('Result.csv',header=True,index=False)


#-----------------------------------##
#群体sample_group 给tree上颜色        ##
##-----------------------------------##
D = {'A':'#E64B35FF','B':'#4DBBD5FF','C':'#00A087FF','D':'#3C5488FF','E':'#F39B7FFF'}
df = pd.read_table(r'E:\中转站\NpX群体\sample_group.txt',sep='\t',names=['sample','group'])
df['color'] = df['group'].map(D)
df['uk'] = 'range'
print(df)
df.iloc[:,[0,3,2,1]].to_csv(r'E:\中转站\NpX群体\sample_group_treecolor.txt',sep='\t',header=False,index=False)


REdf = pd.read_table('prunning.Wholegenome.counts_AAGCTT.txt',sep='\t')
ctgID = []
with open('ZG2ZM_Whole46839ctg.ctgID.txt','r') as f:
    lines = (line.strip() for line in f)
    for line in lines:
        ctgID.append(line)
print(REdf)
print(len(ctgID))
qdf = REdf[REdf['#Contig'].isin(ctgID)]
print(qdf)
qdf.to_csv('ZG2ZM_Whole46839ctg.counts_AAGCTT.txt',sep='\t',header=True,index=False)

import re
import sys

small = []
big = []

with open(sys.argv[1],'r') as f1:
    for line in f1:
        small.append(line.strip())

with open(sys.argv[2],'r') as f2:
    for line in f2:
        big.append(line.strip())

remove = []
for b in big:
    if b not in small:
        remove.append(b)

with open(sys.argv[3],'w') as w:
    for r in remove:
        w.write(r+'\n')

import numpy as np
import matplotlib.pylab as plt

# 构建数据
def model(x, p):
    return x ** (2 * p + 1) / (1 + x ** (2 * p))
x = np.linspace(0.75, 1.25, 201)
# 可视化绘制
fig, ax = plt.subplots(figsize=(4,3),dpi=200)
for p in [10, 15, 20, 30, 50, 100]:
    ax.plot(x, model(x, p), label=p)
    fig.show()



import pandas as pd
import sys
df = pd.read_table(sys.argv[1],sep='\t',names=['firstbin','secondbin','score'])
print(df)
df = df[df['firstbin'] <= 27287]
df = df[df['secondbin'] <=27287]
print(df)
df.to_csv('NpX_100000_iced.modified.noctg.matirx',sep='\t',header=False,index=False)

import sys
from Bio import SeqIO
def parse_fasta(fasta_file):
    D = {}
    for seq in SeqIO.parse(fasta_file,'fasta'):
        D[seq.id] = seq.seq
    N = {}
    for id in D.keys():
        newid = ','.join(id.split('.')[0:2])
        N[newid] = D[id]
    with open('Sb.corrected.cds_primaryTranscriptOnly.fa','w') as w:
        for key in N.keys():
            w.write(key+'\n')
            w.write(str(N[key])+'\n')

if __name__ == '__main__':
    parse_fasta(sys.argv[1])


import pandas as pd
import re
df = pd.read_excel(r'E:\潘浩然\Result\Transit\fpkm_caculate\全套\10.AP85-441.fpkm_dataframe.xlsx')
print(df)
f = lambda x:re.findall('\.([0-9]+)G',x)[0]
df['new'] = df['gene_id'].apply(f)
print(df)


#------------------中果基因组gff3文件geneID的替换----------------------------
import sys
import pandas as pd

def change_name(file):
    D = {
        'Chr01':'Chr1A','Chr02':'Chr1B','Chr03':'Chr1C','Chr04':'Chr1D','Chr05':'Chr1E','Chr06':'Chr1F','Chr07':'Chr1G','Chr08':'Chr1H',
        'Chr09':'Chr1A','Chr10':'Chr1B','Chr11':'Chr1C','Chr12':'Chr1D','Chr13':'Chr1E','Chr14':'Chr1F','Chr15':'Chr1G','Chr16':'Chr1H',
        'Chr17':'Chr1A','Chr18':'Chr18','Chr19':'Chr19','Chr20':'Chr20','Chr21':'Chr21','Chr22':'Chr22','Chr23':'Chr23','Chr24':'Chr1H',
        'Chr25':'Chr1A','Chr26':'Chr1B','Chr27':'Chr1C','Chr28':'Chr1D','Chr29':'Chr1E','Chr30':'Chr1F','Chr31':'Chr31','Chr32':'Chr1H',
        'Chr33':'Chr1A','Chr34':'Chr1B','Chr35':'Chr1C','Chr36':'Chr1D','Chr37':'Chr1E','Chr38':'Chr1F','Chr39':'Chr1G','Chr40':'Chr1H',
        'Chr41':'Chr1A','Chr42':'Chr1B','Chr43':'Chr1C','Chr44':'Chr1D','Chr45':'Chr1E','Chr46':'Chr1F','Chr47':'Chr1G','Chr48':'Chr1H',
        'Chr49':'Chr1A','Chr50':'Chr1B','Chr51':'Chr1C','Chr52':'Chr1D','Chr53':'Chr1E','Chr54':'Chr1F','Chr55':'Chr1G','Chr56':'Chr1H',
        'Chr57':'Chr1A','Chr58':'Chr1B','Chr59':'Chr1C','Chr60':'Chr1D','Chr61':'Chr1E','Chr62':'Chr1F','Chr63':'Chr1G','Chr64':'Chr1H',
        'Chr65':'Chr1A','Chr66':'Chr1B','Chr67':'Chr1C','Chr68':'Chr1D','Chr69':'Chr1E','Chr70':'Chr1F','Chr71':'Chr1G','Chr72':'Chr1H',
        'Chr73':'Chr1A','Chr74':'Chr1B','Chr75':'Chr1C','Chr76':'Chr1D','Chr77':'Chr1E','Chr78':'Chr1F','Chr79':'Chr1G','Chr80':'Chr1H'
         }
    agpdf = pd.read_table(file,names=['seqid','start','end','index','uk1','uk2','uk3','uk4','uk5'])
    agpdf['new'] = agpdf['seqid'].map(D)
    print(agpdf)
    agpdf = agpdf.iloc[:,[-1,1,2,3,4,5,6,7,8]]
    print(agpdf)
if __name__ == '__main__':
    change_name(sys.argv[1])


import pandas as pd
def read_gff(gff):
    D = {
        'Chr01': 'Chr1A', 'Chr02': 'Chr1B', 'Chr03': 'Chr1C', 'Chr04': 'Chr1D', 'Chr05': 'Chr1E', 'Chr06': 'Chr1F',
        'Chr07': 'Chr1G', 'Chr08': 'Chr1H',
        'Chr09': 'Chr1A', 'Chr10': 'Chr1B', 'Chr11': 'Chr1C', 'Chr12': 'Chr1D', 'Chr13': 'Chr1E', 'Chr14': 'Chr1F',
        'Chr15': 'Chr1G', 'Chr16': 'Chr1H',
        'Chr17': 'Chr1A', 'Chr18': 'Chr18', 'Chr19': 'Chr19', 'Chr20': 'Chr20', 'Chr21': 'Chr21', 'Chr22': 'Chr22',
        'Chr23': 'Chr23', 'Chr24': 'Chr1H',
        'Chr25': 'Chr1A', 'Chr26': 'Chr1B', 'Chr27': 'Chr1C', 'Chr28': 'Chr1D', 'Chr29': 'Chr1E', 'Chr30': 'Chr1F',
        'Chr31': 'Chr31', 'Chr32': 'Chr1H',
        'Chr33': 'Chr1A', 'Chr34': 'Chr1B', 'Chr35': 'Chr1C', 'Chr36': 'Chr1D', 'Chr37': 'Chr1E', 'Chr38': 'Chr1F',
        'Chr39': 'Chr1G', 'Chr40': 'Chr1H',
        'Chr41': 'Chr1A', 'Chr42': 'Chr1B', 'Chr43': 'Chr1C', 'Chr44': 'Chr1D', 'Chr45': 'Chr1E', 'Chr46': 'Chr1F',
        'Chr47': 'Chr1G', 'Chr48': 'Chr1H',
        'Chr49': 'Chr1A', 'Chr50': 'Chr1B', 'Chr51': 'Chr1C', 'Chr52': 'Chr1D', 'Chr53': 'Chr1E', 'Chr54': 'Chr1F',
        'Chr55': 'Chr1G', 'Chr56': 'Chr1H',
        'Chr57': 'Chr1A', 'Chr58': 'Chr1B', 'Chr59': 'Chr1C', 'Chr60': 'Chr1D', 'Chr61': 'Chr1E', 'Chr62': 'Chr1F',
        'Chr63': 'Chr1G', 'Chr64': 'Chr1H',
        'Chr65': 'Chr1A', 'Chr66': 'Chr1B', 'Chr67': 'Chr1C', 'Chr68': 'Chr1D', 'Chr69': 'Chr1E', 'Chr70': 'Chr1F',
        'Chr71': 'Chr1G', 'Chr72': 'Chr1H',
        'Chr73': 'Chr1A', 'Chr74': 'Chr1B', 'Chr75': 'Chr1C', 'Chr76': 'Chr1D', 'Chr77': 'Chr1E', 'Chr78': 'Chr1F',
        'Chr79': 'Chr1G', 'Chr80': 'Chr1H'
    }
    colname=['seqid','source','type','start','end','score','strand','phase','attributes']
    df = pd.read_table(gff,comment='#',names=colname)
    df['new'] = df['seqid'].map(D)
    df = df.iloc[:,[-1,1,2,3,4,5,6,7,8]]
    df = df.sort_values(by=['new','start'])
    print(df)
    # prefix = re.findall('([0-9a-zA-Z\_\-\.]+)\.gff3',gff)[0]
    return df

if __name__ == '__main__':
    read_gff(sys.argv[1])


import sys
import pandas as pd

def read_gff(gff):
    colname=['seqid','source','type','start','end','score','strand','phase','attributes']
    df = pd.read_table(gff,comment='#',names=colname)
    # prefix = re.findall('([0-9a-zA-Z\_\-\.]+)\.gff3',gff)[0]

    return df

def read_unanchors(file):
    L = []
    with open(file,'r') as f:
        for line in f:
            L.append(line.strip())
    return L

def gain_unanchorsgff(gffdf,L):
    qdf = gffdf[gffdf['seqid'].isin(L)]
    qdf.to_csv('unanchors.ctg.gff3',sep='\t',header=False,index=False)

if __name__ == '__main__':
    ungffdf = read_gff(sys.argv[1])
    L = read_unanchors(sys.argv[2])
    gain_unanchorsgff(ungffdf,L)


import sys
from Bio import SeqIO

D = {}
with open(sys.argv[1]) as f:
    for line in f:
        line_list = line.strip().split('===>')
        D[line_list[0]] = line_list[1]
# print(D)

D2 = {}
for seq in SeqIO.parse(sys.argv[2], 'fasta'):
    # print(seq.id, seq.seq, sep='\n')
    D2[seq.id] = seq.seq
D3 = {}
for key in D2.keys():
    D3[D[key]] = D2[key]

for key in sorted(D3.items(),key=lambda x:x[1]):
    print(key[0])
    break

import sys
import pandas as pd
import numpy as np
L = []
df = pd.read_excel(sys.argv[1])
df = df.fillna('0')
for index,row in df.iloc[:,3:].iterrows():
    #print(row)
    for i in row:
        #print(i)
        if i.startswith('SoZg') and ',' not in i:
            L.append(i)
        elif i.startswith('SoZg') and ',' in i:
            line_list = i.split(',')
            for x in line_list:
                L.append(x)
with open(sys.argv[2],'w') as w:
    for l in L:
        w.write(l+'\n')

