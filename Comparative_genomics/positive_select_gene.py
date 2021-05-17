#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/5/14


import re
from path import Path
import argparse
import sys
sys.path.append('/public1/home/stu_panhaoran/biosoft/pyhr/scripts')
import runshell
from rich.progress import track


__author__ = 'Haoran Pan'
__mail__ = 'haoran_pan@qq.com'
__date__ = '20210516'
__version__ = 'v1.0'


def filter_cds(phylip_Path,output_Path):
    num =0
    files = [file for file in Path(phylip_Path).files() if file.endswith('.phylip')]
    for file in track(files):
        file_name = re.findall('(OG[0-9]+)\.phylip',file.basename())[0]

        with open(file,'r') as f:
            f = f.readline()
            specie_number, seq_len = f.split('\t')
            if int(seq_len) % 3 == 0:
                num +=1
                null_read_content = []
                runshell.run_command_noprint('cp OGnumber.null.ctl {}.null.read.ctl'.format(file_name))
                with open('{}.null.read.ctl'.format(file_name),'r') as f1:
                    # lines = (line.strip() for line in f1)
                    for line in f1:
                        null_read_content.append(line.rstrip())
                with open('{}.null.ctl'.format(file_name),'w') as f2:
                    null_write_content = []
                    for row in null_read_content:
                        if 'seqfile' in row:
                            row = re.sub('path',phylip_Path,row)
                            row = re.sub('OGnumber',file_name,row)
                            null_write_content.append(row)
                        elif 'outfile' in row:
                            row = re.sub('path', output_Path, row)
                            row = re.sub('OGnumber',file_name,row)
                            null_write_content.append(row)
                        else:
                            null_write_content.append(row)
                    f2.write('\n'.join(null_write_content))

                alte_read_content = []
                runshell.run_command_noprint('cp OGnumber.alte.ctl {}.alte.read.ctl'.format(file_name))
                with open('{}.alte.read.ctl'.format(file_name),'r') as f3:
                    # lines = (line.strip() for line in f1)
                    for line in f3:
                        alte_read_content.append(line.rstrip())
                with open('{}.alte.ctl'.format(file_name),'w') as f2:
                    alte_write_content = []
                    for row in alte_read_content:
                        if 'seqfile' in row:
                            row = re.sub('path',phylip_Path,row)
                            row = re.sub('OGnumber',file_name,row)
                            alte_write_content.append(row)
                        elif 'outfile' in row:
                            row = re.sub('path', output_Path, row)
                            row = re.sub('OGnumber',file_name,row)
                            alte_write_content.append(row)
                        else:
                            alte_write_content.append(row)
                    f2.write('\n'.join(alte_write_content))
                runshell.run_command_noprint('rm {}.alte.read.ctl {}.null.read.ctl'.format(file_name,file_name))
    runshell.run_command('rm OGnumber.null.ctl OGnumber.alte.ctl')
    print(num)


def store_lnL(result_path):
    alte_files = [file for file in Path(result_path).files() if file.endswith('.alte.txt')]
    null_files = [file for file in Path(result_path).files() if file.endswith('.null.txt')]
    alte2lnL = {}
    for alte in alte_files:
        res = runshell.run_command_captureoutput('grep "lnL" {}'.format(alte))
        if res:
            alte2lnL[re.findall('(OG[0-9]+)\.alte',alte)[0]] = float(re.findall('[\-\.0-9]+',res)[2])

    null2lnL = {}
    for null in null_files:
        res = runshell.run_command_captureoutput('grep "lnL" {}'.format(null))
        if res:
            null2lnL[re.findall('(OG[0-9]+)\.null',null)[0]] = float(re.findall('[\-\.0-9]+',res)[2])
    return alte2lnL,null2lnL


def cal_2l(alte2lnL_Dict,null2lnL_Dict):
    OGnumber22l = {}
    for key in alte2lnL_Dict.keys():
        o2l = 2 * abs(alte2lnL_Dict[key] - null2lnL_Dict[key])
        OGnumber22l[key] = o2l
    return OGnumber22l


def chi_test(OGnumber22l_Dict,result_Path):
    with open('{}//PSG_result.tab'.format(result_Path),'w') as f:
        for key in OGnumber22l_Dict.keys():
            res = runshell.run_command_captureoutput('/public1/home/stu_panhaoran/biosoft/paml/paml4.9j/src/chi2 1.5 {}'.format(OGnumber22l_Dict[key]))
            chi_value = float([tmp for tmp in res.split('\n')][1].split('=')[2])/2
            if chi_value <= 0.05 :
                f.write('{}{}{}{}{}{}'.format(key,'\t',chi_value,'\t','Yes','\n'))
            elif chi_value > 0.05 :
                f.write('{}{}{}{}{}{}'.format(key,'\t',chi_value,'\t','No','\n'))


def main(args):
    cds_msa_phy_path = args.phypath
    ctl_output_path = args.ctloutpath
    mode = args.mode
    if mode == 'filter':
        filter_cds(cds_msa_phy_path,ctl_output_path)
    elif mode == 'cal':
        alte_dict, null_dict = store_lnL(ctl_output_path)
        OGnumber22l = cal_2l(alte_dict,null_dict)
        chi_test(OGnumber22l, ctl_output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''This script uses paml to filter positively selected genes''',
        usage="python {} -m filter --phypath --ctloutpath \n or \npython {} -m cal --ctloutpath ".format(sys.argv[0],sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__, __mail__, __date__,
                                                                            __version__))
    parser.add_argument('--phypath', required=False, help='Input cds_msa_phylip path')
    parser.add_argument('--ctloutpath', required=False, help='Specify the ctl.txt file output path')
    parser.add_argument('-m','--mode', required=True,choices = ['filter','cal'], help='choose run mode')
    args = parser.parse_args()
    main(args)