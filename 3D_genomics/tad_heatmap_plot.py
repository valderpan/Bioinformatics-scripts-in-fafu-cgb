# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-


import subprocess
import datetime
import argparse

chromsize_dict = {}
def get_chrom_sizes(chromsizes_file):
    with open(chromsizes_file) as f:
        for line in f:
            line = line.strip()
            chrom,size = line.split()
            chromsize_dict[chrom] = size

    return chromsize_dict

def run_command(cmd):
    print(cmd)
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
        print("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd,datetime.datetime.now()))

def run_tad_plot(speciesname,cool,hitad_out):
    for i in chromsize_dict.keys():
        cmd = 'tad-plot -O '+speciesname+'.'+i+'.pdf -p '+cool+' -T '+hitad_out+' -C '+i+' -W RAW -S 0 -E '+chromsize_dict[i]
        run_command(cmd)

def main(args):
    chromsizes_file = args.chromsize
    speciesname = args.speciesname
    cool =args.cool
    hitad_out = args.hitadout
    get_chrom_sizes(chromsizes_file)
    run_tad_plot(speciesname,cool,hitad_out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='tad_heatmap_plot.py',
        description='''Draw the heat map of tad at the chromosome level according to the tad-plot script in hitad''')
    parser.add_argument('-cs', '--chromsize', required=True, help=' Input chromsizes file ')
    parser.add_argument('-sn', '--speciesname', required=True, help='Input species name')
    parser.add_argument('-co', '--cool', required=True, help='Input cool file')
    parser.add_argument('-hi', '--hitadout', required=True, help='Input hitad output file')
    args = parser.parse_args()
    main(args)


