# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

import pandas as pd


def reads_mapping_calcultate(Ks_file1,ks_name,Upper_limit):
    df1 = pd.read_table(Ks_file1,sep='\t')
    ks_df1 = df1.iloc[:,[0,3]]
    ks_df1 = ks_df1[(ks_df1['dS-ng']>=0) & (ks_df1['dS-ng']<Upper_limit)]
    ks_df1.iloc[:,1] = round(ks_df1.iloc[:,1],2)
    ks_count_dict1 = ks_df1.iloc[:,1].value_counts().to_dict()
    ks_count_df1 = pd.DataFrame(ks_count_dict1,index = [''.join(ks_name)+'_counts']).T
    ks_count_df1[''.join(ks_name)] = ks_count_df1[''.join(ks_name)+'_counts']/sum(ks_count_df1[''.join(ks_name)+'_counts'])*100
    output = ks_count_df1.sort_index()
    output = output.reset_index()
    return  output

