import os
import pandas as pd


os.chdir(r'/public1/home/stu_panhaoran/CRE_biotask/v0806/split_1M/plantCARE_output/onlygroupby')
sub_file_list = os.listdir()
sub_file = []
for i in sub_file_list:
    # print(i)
    df = pd.read_table(i,sep='\t',header=0)
    # print(df)
    sub_file.append(df)
    
df2 = pd.concat(sub_file)
#df2.to_csv('/public1/home/stu_panhaoran/CRE_biotask/concat_file.csv',index=False)
df2.to_csv('/public1/home/stu_panhaoran/CRE_biotask/v0806/Sspon_all_genome_CRE.tab',sep='\t',header=True,index=False)
