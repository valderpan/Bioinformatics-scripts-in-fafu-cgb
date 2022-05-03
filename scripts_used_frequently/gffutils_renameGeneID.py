#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/05/03

import gffutils
import argparse
import sys
import os.path as op
from collections import OrderedDict


def createDB(gff, dbfn=None):
    """
    create a db file for gffutils
    Params:
    --------
    gff: `str` gff file
    dbfn: `str` db name [default: gff.db]
    Returns:
    --------
    db: `gffutils.FeatureDB` db for gffutils
    Examples:
    ---------
    >>> createDB('in.gff3')
    <FeatureDB ...>
    """
    ## create gff db
    if not dbfn: 
        db_name = gff + ".db"
    else:
        db_name = dbfn
    if not op.exists(db_name):
        print('No such database file of `{}`, creating ...'.format(db_name))
        gffutils.create_db(gff, dbfn=db_name, keep_order=True)
    else:
        print('Already exists DB file of `{}`, skip.'.format(db_name))
    db = gffutils.FeatureDB(db_name)
    
    return db
    
def main(args):
    output = args.output

    db = createDB(args.gff, args.dbname)
    rename_db = OrderedDict(i.strip().split() for i in open(args.renameList)
                        if i.strip())
    
    for ID in rename_db:
        print(str(db[ID]).replace(ID, rename_db[ID]), file=output)
        for feature in db.children(ID, order_by='start'):
            print(str(feature).replace(ID, rename_db[ID]), file=output)
        print("", file=output)
    print("Successful. Output file is in `{}`".format(output.name))



if __name__ == "__main__":
    p = argparse.ArgumentParser(prog=sys.argv[0],
                        description='''This script is used to rename ID in gff file (given gene ID, it can also automatically replace the ID of its CDS, exon, etc.)''',
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('gff', 
            help='gff file ')
    pReq.add_argument('renameList', 
            help='rename list, two columns <old_gene_name\tnew_gene_name>')
    pOpt.add_argument('-d', '--dbname', default=None,
            help='gff database name [default: gff3.db]')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    args = p.parse_args()
    main(args)