#!/usr/bin/env python

# for tgi cluster:
#/gapp/x64linux/opt/pythonbrew/venvs/Python-2.7.6/gemini/bin/python
# for uva cluster:

import pysam
import sys
import argparse
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-01-01 16:58 $"

def bamfilterrg(bamfile, readgroup, limit, is_sam, bam_out, uncompressed_out):
    # set input file
    if bamfile is None: 
        if is_sam:
            in_bam = pysam.AlignmentFile('-', 'r')
        else:
            in_bam = pysam.AlignmentFile('-', 'rb')
    else:
        if is_sam:
            in_bam = pysam.AlignmentFile(bamfile, 'r')
        else:
            in_bam = pysam.AlignmentFile(bamfile, 'rb')

    # set output file
    if uncompressed_out:
        out_bam = pysam.AlignmentFile('-', 'wbu', template=in_bam)
    elif bam_out:
        out_bam = pysam.AlignmentFile('-', 'wb', template=in_bam)
    else:
        out_bam = pysam.AlignmentFile('-', 'wh', template=in_bam)
        

    # parse readgroup string
    try:
        rg_list = readgroup.split(',')
    except AttributeError:
        rg_list = None

    counter = 0
    for al in in_bam:
        # must be in a user specified readgroup
        if rg_list and al.get_tag('RG') not in rg_list:
            continue

        # write out alignment
        out_bam.write(al)
        counter += 1
        
        # bail if reached limit
        if limit is not None and counter >= limit:
            break

# ============================================
# functions
# ============================================

# class that holds reads from a sequence fragment
class Namegroup():
    def __init__(self, al):
        self.alignments = list()
        self.name = al.query_name
        self.sa = 0
        self.num_prim = 0
        self.add_alignment(al)

    def add_alignment(self, al):
        self.alignments.append(al)
        if not al.is_secondary:
            self.num_prim += 1
            try:
                # print self.sa
                self.sa += len(al.get_tag('SA').rstrip(';').split(';'))
            except KeyError:
                pass

    def is_complete(self):
        return self.num_prim == 2 and len(self.alignments) == self.sa + 2

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
bamfilterrg.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: filter readgroup(s) from a BAM file")
    parser.add_argument('-i', '--input', metavar='BAM', required=False, help='Input BAM file')
    parser.add_argument('-r', '--readgroup', metavar='STR', default=None, required=False, help='Read group(s) to extract (comma separated)')
    parser.add_argument('-n', metavar='INT', type=int, default=None, required=False, help='Output first n alignments and quit')
    parser.add_argument('-S', required=False, action='store_true', help='Input is SAM format')
    parser.add_argument('-b', required=False, action='store_true', help='Output BAM format')
    parser.add_argument('-u', required=False, action='store_true', help='Output uncompressed BAM format (implies -b)')

    # parse the arguments
    args = parser.parse_args()

    # bail if no BAM file
    if args.input is None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
    
    # send back the user input
    return args

# ============================================
# driver
# ============================================

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main():
    args = get_args()
    bamfilterrg(args.input, args.readgroup, args.n, args.S, args.b, args.u)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError as e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
    
