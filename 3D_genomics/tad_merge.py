# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-


import argparse
from intervaltree import IntervalTree, Interval

def get_smaller_tad(hitad_outfile,hitad_merged_outfile):
    tree_dict = {}
    with open(hitad_outfile) as f1:
        for line in f1 :
            line = line.strip()
            chrom,start,end,rank = line.split()
            start,end = int(start),int(end)
            if chrom not in tree_dict:
                tree_dict[chrom] = IntervalTree()
                tree_dict[chrom].addi(start, end, rank)
            overlap = list(tree_dict[chrom].overlap(start, end))
            if overlap:
                if (end - start) < overlap[0].length():
                    tree_dict[chrom].remove(overlap[0])
                    tree_dict[chrom].addi(start, end, rank)
            else:
                tree_dict[chrom].addi(start, end, rank)
    with open(hitad_merged_outfile, 'w') as f2:
        for chrom in tree_dict.keys():
            for item in sorted(tree_dict[chrom]):
                f2.write("\t".join((chrom, str(item.begin), str(item.end), item.data))+'\n')
def main(args):
    hitad_outfile = args.input
    hitad_merged_outfile = args.output
    get_smaller_tad(hitad_outfile,hitad_merged_outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='tad_merge.py',
        description='''Merge TAD to find smaller TAD''')
    parser.add_argument('-i','--input',required=True,help='Input hiad output file')
    parser.add_argument('-o','--output',required=True,help='Output the merged file')
    args = parser.parse_args()
    main(args)
