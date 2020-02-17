# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

# 合并目录下的所有excel表格
import os
import pandas as pd

os.chdir(r'D:\Result\Transit\fpkm_caculate\单套\ses昼夜节律')

file_list = os.listdir('./')
# print(file_list)

first_file = pd.read_excel('1.xlsx')
# print(first_file)
file_dict = {}
for i in file_list[1:]:
    file_dict[i] = pd.read_excel(i)
    first_file = pd.concat([first_file,file_dict[i]])
first_file.to_excel('ses叶段.xlsx',index=False,header=True)