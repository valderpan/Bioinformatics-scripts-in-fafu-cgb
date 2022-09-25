#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/09/23

'''
%prog <gff3 file> <ctgID file> <gene bed file(jcvi output)> <output file>

Copy the corresponding contig gff infomation
After collapse_rescue, a new contig with the _d suffix will be generated, 
and the new contig will be assigned a value according to the previous contig gff information

>>> python %prog gff3 ctgID.txt gene.bed output

__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20220924'
__version__ = 'v1.0'

'''


import re
import sys
import pandas as pd
from rich.traceback import install


def import_gff(gff_path, _type=None):
    compression = 'gzip' if gff_path.endswith('.gz') else 'infer'
    gff_df = pd.read_csv(gff_path, sep='\t', 
                         header=None, 
                         comment="#", 
                         compression = compression,
                         names=('seqid', 'source',
                                'type', 'start', 
                                'end', 'score', 
                                'strand','phase',
                                'attributes'))
    if _type: 
        if _type in set(gff_df.head(100)['type']):
            gff_df = gff_df[gff_df['type'] == _type]
        else:
            print('Warning: Failed to filter data by {}, input type is not correct'.format(_type))
    return gff_df 

def read_ctg(ctgId):
    ctgL = []
    with open(ctgId) as f:
        lines = (line.strip() for line in f)
        for line in lines:
            ctgL.append(line)
    return ctgL

def read_bed(bed):
    beddf = pd.read_table(bed,sep='\t',names=['seqid','start','end','geneid','uk1','uk2'])
    seq2gene = {}
    for i in beddf.groupby('seqid'):
        seq2gene[i[0]] = i[1]['geneid'].tolist()
    return seq2gene

def concatLL(gffdf,ctgL,seq2gene,output):
    concatL = []
    gffctg = gffdf['seqid'].unique().tolist()
    for index,i in enumerate(ctgL):
        if i in gffctg:
            qdf = gffdf[gffdf['seqid']==i]
            concatL.append(qdf)
        else:
            if '_d' in i:
                newi = i.split('_')[0]
                qdf = gffdf[gffdf['seqid']==newi]
                # print(qdf)
                # print('='*20)
                # qdf_gene = qdf[qdf['type'=='gene']]
                # print(newi)
                if newi in seq2gene.keys():
                    for j in seq2gene[newi]:
                        newgeneid = int(re.findall('ZG([0-9]+)',j)[0])+index+247678
                        # if len(newgeneid) < 6:
                        #     newgeneid = '0{}'.format(newgeneid)
                        qdf['attributes'] = qdf['attributes'].str.replace(j,'ZG{}'.format(newgeneid))
                    qdf['seqid'] = qdf['seqid'].str.replace(newi,i)
                    # print(qdf)
                    concatL.append(qdf)
                else:
                    print('{} to {} not in ctg.gff3'.format(i,newi))
    res = pd.concat(concatL)
    res.to_csv(output,sep='\t',header=False,index=False)

# if __name__ == '__main__':
#     install()
#     gffdf= import_gff(sys.argv[1])
#     # print(gffdf)
#     ctgL = read_ctg(sys.argv[2])
#     seq2gene = read_bed(sys.argv[3])
#     concatLL(gffdf,ctgL,seq2gene,sys.argv[4])

if __name__ == "__main__":
    from optparse import OptionParser
    from rich.traceback import install
    install()
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args)==4:
        gff_file = args[0]
        ctg_file = args[1]
        bed_file = args[2]
        output = args[3]
        gffdf= import_gff(gff_file)
        ctgL = read_ctg(ctg_file)
        seq2gene = read_bed(bed_file)
        concatLL(gffdf,ctgL,seq2gene,output)
    else:
        sys.exit(p.print_help())