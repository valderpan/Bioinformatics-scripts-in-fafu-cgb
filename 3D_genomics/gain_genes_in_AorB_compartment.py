# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-


'''
需求：在喜马拉雅割手密每条染色体上找A/B compartment各自区域内的基因，查看其表达量情况
'''

import pandas as pd
import os
import re
from interval import Interval
import math
import argparse


def get_chr_AB_region(Hitc_result_path,output_path1):
    '''
    此函数根据R包HiTC的输出结果文件，将整条染色体划分为A/B区域，并将每条染色体上的分别输出
    input_path = r'D:\\Result\\TAD\\Nepal_AP85_compartment_0715'
    output_path = r'C:\\Users\\dell\\Desktop\\Nepal_project\\Nepal_chr_AB_geneexpression'
    '''
    path1 = Hitc_result_path
    os.chdir(path1)
    path1_list = os.listdir(path1)
    for i in path1_list:
        compartment_df = pd.read_csv(i,header=0)
        compartment_df = compartment_df.iloc[:,[2,3,6]]
        concat_A_list = []
        concat_B_list = []
        for index,row in compartment_df.iterrows():
            if row['score'] > 0 :
                concat_A_list.append(compartment_df.loc[[index],:])
            else:
                concat_B_list.append(compartment_df.loc[[index], :])
        A_compament_df = pd.concat(concat_A_list)
        B_compament_df = pd.concat(concat_B_list)
        A_region = A_compament_df.iloc[:,0:2]
        output_name = re.findall('([A-Z0-9a-z\_]+)\.',i)
        A_region.to_csv(output_path1+'\\'+output_name[0]+'_A.region.txt',sep='\t',header=True,index=False)
        B_region = B_compament_df.iloc[:, 0:2]
        B_region.to_csv(output_path1 + '\\' + output_name[0] + '_B.region.txt', sep='\t',header=True,index=False)


def get_gene_position(gff3_file,output_path2,skiprow=1):
    '''
    此函数对gff3文件进行处理，去掉无用信息，将Chr,gene,start,end,ID 5列进行输出
    gff3_file:r'D:\\Result\\Nipor_Ss.v20200409.Chr.gff3'
    output_path2:r'C:\\Users\\dell\\Desktop\\Nepal_project\\Nepal_chr_AB_geneexpression'
    '''
    gff_df = pd.read_table(gff3_file,sep='\t',skiprows=skiprow,header=None)
    gene_df = gff_df[(gff_df.iloc[:,2].str.contains('gene')) & (gff_df.iloc[:,0].str.contains('Chr'))] #找到gene的注释行并将contig的注释行删除
    gene_df = gene_df.iloc[:,[0,2,3,4,8]]    #提取有用信息列(chr,gene,start,end,Attributes)
    Attributes = gene_df.iloc[:,4].str.split(';',expand=True)
    Attributes = Attributes.iloc[:,0].str.split('=',expand=True)
    gene_df = pd.merge(gene_df,Attributes.iloc[:,1],on=gene_df.index)
    gene_df = gene_df.iloc[:,[1,2,3,4,6]]      #提取merge后的有用信息列(Chr,gene,start,end,ID)
    gene_df.to_csv(output_path2 +'\\'+'geneID_position.txt',sep='\t',header=None,index=False)


def match_ABregion2gene_position(AB_region_file_path,gff3_file,output_path3):
    '''
    此函数是将上一步生成的A/B compartment region文件，根据处理后的gff文件上每个基因的位置，将在region上的geneID输出
    考虑到有些基因一半落在A compartment区域，一半落在B compartment区域，这里根据在两部分的占比大小来决定基因最终的归属
    这一步最终生成A/B compartment区域上有哪些基因的文件
    AB_region_file_path:r'C:\\Users\\dell\\Desktop\\Nepal_project\\Nepal_chr_AB_geneexpression'
    gff3_file:r'C:\\Users\\dell\\Desktop\\Nepal_project\\Nepal_chr_AB_geneexpression\\geneID_position.txt'
    output_path3:r'C:\\Users\\dell\\Desktop\\Nepal_project\\Nepal_chr_AB_geneexpression\\Total_gene_filted'
    '''
    os.chdir(AB_region_file_path)
    file_list = os.listdir(AB_region_file_path)
    # print(file_list)
    for file in file_list:           #这里的每个file都是500kb的bins
        outputfile_want_name = re.findall('([A-Z0-9a-z\_]+)\.',file)[0]
        print(outputfile_want_name)
        file_chr_name = re.findall('([A-Za-z0-9]+)\_',file)[1]
        # print(file_chr_name)
        AorB_name  = re.findall('([A-Z]+)\.',file)[0]
        # print(AorB_name)
        AB_region_inter_list = []
        with open(file) as f:
            for line in f :
                line = line.strip()
                if line.startswith('start'):
                    continue
                inter_start,inter_end = line.split()
                inter = Interval(int(inter_start),int(inter_end))
                AB_region_inter_list.append(inter)
        gene_df = pd.read_table(gff3_file, sep='\t', header=None)
        gene_df.rename(columns = {0:'Chr',1:'gene',2:'start',3:'end',4:'ID'},inplace=True)
        gene_df = gene_df[gene_df['Chr']==file_chr_name]
        if AorB_name == 'A':
            A_gene_list = []
            A_total_list = []
            result = open(output_path3+'\\'+outputfile_want_name+'.genes.txt','w')
            for index,row in gene_df.iterrows():
                gene_inter_start = int(gene_df.loc[index]['start'])
                gene_inter_end = int(gene_df.loc[index]['end'])
                gene_inter = Interval(gene_inter_start,gene_inter_end)
                for i in AB_region_inter_list:
                    if gene_inter.overlaps(i):
                        A_total_list.append(gene_df.loc[index]['ID'])
                        if gene_inter.lower_bound <= i.lower_bound <= gene_inter.upper_bound:
                            if (i.lower_bound-gene_inter.lower_bound) < (gene_inter.upper_bound-i.lower_bound):
                                A_gene_list.append(gene_df.loc[index]['ID'])
                        elif gene_inter.lower_bound <= i.upper_bound <= gene_inter.upper_bound:
                            if (i.upper_bound-gene_inter.lower_bound) > (gene_inter.upper_bound-i.upper_bound):
                                A_gene_list.append(gene_df.loc[index]['ID'])
                        elif gene_inter.lower_bound >= i.lower_bound and gene_inter.upper_bound <= i.upper_bound:
                            A_gene_list.append(gene_df.loc[index]['ID'])
            print('A compartment filted gene number is {}'.format(len(A_gene_list)))
            for gene in A_gene_list:
                result.write(gene+'\n')
            result.close()
            print('A compartment total gene number(nofilted) is {}'.format(len(A_total_list)))

        elif AorB_name == 'B':
            B_gene_list = []
            B_total_list = []
            result = open(output_path3 + '\\' + outputfile_want_name + '.genes.txt', 'w')
            for index,row in gene_df.iterrows():
                gene_inter_start = int(gene_df.loc[index]['start'])
                gene_inter_end = int(gene_df.loc[index]['end'])
                gene_inter = Interval(gene_inter_start,gene_inter_end)
                for i in AB_region_inter_list:
                    if gene_inter.overlaps(i):
                        B_total_list.append(gene_df.loc[index]['ID'])
                        if gene_inter.lower_bound <= i.lower_bound <= gene_inter.upper_bound:
                            if (i.lower_bound-gene_inter.lower_bound) < (gene_inter.upper_bound-i.lower_bound):
                                B_gene_list.append(gene_df.loc[index]['ID'])
                        elif gene_inter.lower_bound <= i.upper_bound <= gene_inter.upper_bound:
                            if (i.upper_bound-gene_inter.lower_bound) > (gene_inter.upper_bound-i.upper_bound):
                                B_gene_list.append(gene_df.loc[index]['ID'])
                        elif gene_inter.lower_bound >= i.lower_bound and gene_inter.upper_bound <= i.upper_bound:
                            B_gene_list.append(gene_df.loc[index]['ID'])
            print('B compartment filted gene number is {}'.format(len(B_gene_list)))
            for gene in B_gene_list:
                result.write(gene+'\n')
            result.close()
            print('B compartment total gene number(nofilted) is {}'.format(len(B_total_list)))



def extract_gene_expression_value(fpkm_file,gene_file_path,output_path4):
    '''
    此函数用于将A/B compartment 区域上存在的基因的表达量的统计.
    根据fpkm表达谱，将各组织表达量取log2并计算表达值。(这里取得是L叶片和S茎)
    fpkm_file:r'C:\\Users\\dell\\Desktop\\Nepal_project\\Nepal_fpkm_df.raw.xlsx'
    gene_file_path:r'C:\\Users\\dell\\Desktop\\Nepal_project\\Nepal_chr_AB_geneexpression\\Total_gene_filted'
    output_path4:r'C:\\Users\\dell\\Desktop\\Nepal_project\\Nepal_chr_AB_geneexpression\\Total_gene_filted\\Total_gene_fpkm'
    '''
    fpkm_df = pd.read_excel(fpkm_file)
    fpkm_df['L'] = (fpkm_df.iloc[:,1]+fpkm_df.iloc[:,2]+fpkm_df.iloc[:,3])/3
    fpkm_df['S'] = (fpkm_df.iloc[:,4]+fpkm_df.iloc[:,5]+fpkm_df.iloc[:,6])/3
    fpkm_df = fpkm_df.iloc[:,[0,-2,-1]]
    os.chdir(gene_file_path)
    file_list = os.listdir(gene_file_path)
    # print(file_list)
    average_L_dict = {}
    average_S_dict = {}
    A_gene_num = 0
    B_gene_num = 0
    for file in file_list[:-1]:
        print(file)
        file_chr_name = re.findall('([A-Z0-9a-z\_]+)\.', file)[0]
        print(file_chr_name)
        with open(file) as f:
            gene_list = [line.strip() for line in f]
            if file_chr_name[-1] == 'A':
                A_gene_num += len(gene_list)
            elif file_chr_name[-1] == 'B':
                B_gene_num += len(gene_list)
            concat_list = []
            for index,row in fpkm_df.iterrows():
                if fpkm_df.loc[index]['gene_id'] in gene_list:
                    concat_list.append(fpkm_df.loc[[index]])
            result = pd.concat(concat_list)
            result['L1'] = result['L'].apply(lambda x:math.log2(x+1))
            result['S1'] = result['S'].apply(lambda x:math.log2(x+1))
            result = result.iloc[:,[0,3,4]]
            # print(result)
            L_average = result['L1'].mean()
            S_average = result['S1'].mean()
            if file_chr_name[-1] == 'A':
                average_L_dict[file_chr_name] = round(L_average,2)
                average_S_dict[file_chr_name] = round(S_average,2)
            elif file_chr_name[-1] == 'B':
                average_L_dict[file_chr_name] = round(L_average,2)
                average_S_dict[file_chr_name] = round(S_average,2)

    L_average_fpkm = pd.DataFrame(average_L_dict,index=[0])
    S_average_fpkm = pd.DataFrame(average_S_dict,index=[0])
    print(L_average_fpkm.T)
    print(S_average_fpkm.T)
    L_average_fpkm.T.to_csv(output_path4+'\\'+'L_fpkm.txt',sep='\t',header=True,index=True)
    S_average_fpkm.T.to_csv(output_path4+'\\'+'S_fpkm.txt',sep='\t',header=True,index=True)
    print('Total A compartment region gene number is {}'.format(A_gene_num))
    print('Total B compartment region gene number is {}'.format(B_gene_num))


def main(args):
    if args.hitcresultpath and args.ABregionpath and args.gff and args.gffresultpath:
        Hitc_result_path = args.hitcresultpath
        output_path1 = args.ABregionpath
        gff3_file = args.gff
        output_path2 = args.gffresultpath
        skiprow = args.skiprow
        get_chr_AB_region(Hitc_result_path,output_path1)
        get_gene_position(gff3_file,output_path2,skiprow)
    if args.geneposition and args.geneIDpath and args.ABregionpath:
        ABregion_path = args.ABregionpath
        gene_position = args.geneposition
        output_path3 = args.geneIDpath
        match_ABregion2gene_position(ABregion_path,gene_position,output_path3)
    if args.fpkmfile and args.geneIDpath and args.fpkmresultpath:
        fpkm_file = args.fpkmfile
        gene_file_path = args.geneIDpath
        output_path4 = args.fpkmresultpath
        extract_gene_expression_value(fpkm_file,gene_file_path,output_path4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='gain_genes_in_AorB_compartment.py',
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Find out which genes are in the A/B compartment region and check the expression level of the genes
        step 1
        python gain_genes_in_AorB_compartment.py -hp HitcResultPath -rg ABRegionPath -g gff -gp GffResultPath
        
        step 2
        python gain_genes_in_AorB_compartment.py -rg GffResultPath -ge GenePosition -gi GeneIDPath
        
        step 3
        python gain_genes_in_AorB_compartment.py -ff FpkmFile -gi GeneIDPath -fp FpkmResultPath
        ''')
    parser.add_argument('-hp', '--hitcresultpath', help="Input Hitc's compartment statistics result path")
    parser.add_argument('-rg', '--ABregionpath', help="The bins file obtained according to the result of Hitc's compartment")
    parser.add_argument('-g', '--gff',  help='Input gff annotation file')
    parser.add_argument('-sk', '--skiprow',type=int  ,help='The header that needs to be skipped in the gff file ')
    parser.add_argument('-gp', '--gffresultpath',  help='Output gene location file path')
    parser.add_argument('-ge', '--geneposition',  help='Input gene location file')
    parser.add_argument('-gi', '--geneIDpath',  help='Output the genes in the A/B region on each chromosome')
    parser.add_argument('-ff', '--fpkmfile',  help='Input RNA-Seq gene expression file')
    parser.add_argument('-fp', '--fpkmresultpath',  help='Output genes average fpkm value in A/B region on each chromosome')
    args = parser.parse_args()
    main(args)

