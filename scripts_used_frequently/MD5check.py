#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/5/11

'''
%prog [moved file path] <origin_md5_file>

Use the md5 value of the file in the moved directory to compare with the md5 of the file stored in the original directory
'''


import sys
import subprocess
sys.path.append('/public1/home/stu_panhaoran')
import SetLog
import re


def cal_MD5(file_path):
    movedfile2md5 = {}
    if file_path[-1] == '/':
        res = subprocess.check_output('md5sum {}*.gz'.format(file_path),shell=True)
        res = res.decode('utf-8')

        res_list = [tmp.split('  ') for tmp in res.split('\n')]

        for i in res_list:
            if len(i[0]) > 1:
                name = [i for i in re.findall('\/([0-9a-zA-Z\-\_\.]+)', i[1]) if i.endswith('.gz')][0]
                movedfile2md5[name] = i[0]
    else:
        res = subprocess.check_output('md5sum {}/*.gz'.format(file_path), shell=True)
        res = res.decode('utf-8')
        res_list = [tmp.split('  ') for tmp in res.split('\n')]
        for i in res_list:
            if len(i[0]) > 1:
                name = [i for i in re.findall('\/([0-9a-zA-Z\-\_\.]+)',i[1]) if i.endswith('.gz')][0]
                movedfile2md5[name] = i[0]
    return movedfile2md5


def Origin_md5(origin_md5_file):
    file2md5 = {}
    with open(origin_md5_file,'r') as f:
        lines = (line.strip() for line in f)
        for line in lines:
            md5,name = line.split('  ')
            name = [i for i in re.findall('\/([0-9a-zA-Z\-\_\.]+)',name) if i.endswith('.gz')][0]
            file2md5[name] = md5
    return file2md5


def check_MD5(Originfile2md5,Movefile2md5):
    if len(Originfile2md5) != len(Movefile2md5):
        SetLog.error_out('Inconsistent number of files !')
    else:
        SetLog.info_out('File number is OK ~')
        for ori in Originfile2md5.keys():
            if Movefile2md5[ori] == Originfile2md5[ori]:
                SetLog.info_out('File {} MD5 is OK ~'.format(ori))
            else:
                SetLog.error_out('File {} MD5 value is wrong !!'.format(ori))



if __name__ == "__main__":
    from optparse import OptionParser
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args)==2:
        moved_file_path = args[0]
        origin_md5_file = args[1]
        movedfile2md5 = cal_MD5(moved_file_path)
        orifinfile2md5 = Origin_md5(origin_md5_file)
        check_MD5(movedfile2md5,orifinfile2md5)
    else:
        sys.exit(p.print_help())
