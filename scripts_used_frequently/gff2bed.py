#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/4/8

import sys
import argparse

__author__ = 'Haoran Pan'
__mail__ = 'haoran_pan@qq.com'
__date__ = '20210408'
__version__ = 'v1.0'

def gff2bed(gfffile,feature='gene',target='ID'):
    with open(gfffile,'r') as f:
        lines = (line.strip() for line in f)
        for line in lines:
            if line.startswith('#'):
                continue
            else:
                line_list = line.split('\t')
                if line_list[2] == feature:
                    geneID = line_list[8].split(';')[0].split(target+'=')[1]
                    out = [line_list[0],line_list[3],line_list[4],geneID]
                    print("\t".join(map(str, out)))

def main(args):
    gff = args.gfffile
    feature = args.feature
    target = args.target
    gff2bed(gff,feature,target)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Convert GFF files to bed files''',
        usage="python {} -g gff -f feature -t target > .bed".format(sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__, __mail__, __date__,__version__))
    parser.add_argument('-g', '--gfffile', dest='gfffile',required=True, help='Input gff file')
    parser.add_argument('-f', '--feature', dest='feature',default='gene',required=False, help='extract feature of annotations file [default: gene]')
    parser.add_argument('-t', '--target', dest='target',default='ID',required=False, help='the target of feature [default: ID]')
    args = parser.parse_args()
    main(args)
