#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on  2022/4/23 22:35

@author: quinn
"""
import pysam
import sys
from collections import defaultdict

sam_file = sys.argv[1]
cell_filter = sys.argv[2]
out_file = sam_file.split('.')[0] + '_barcode_and_gene.txt'

cell = defaultdict()
with open(cell_filter, 'r') as f:
    for line in f:
        line = line.strip()
        cell[line] = '0'

with pysam.AlignmentFile(sam_file, 'r') as f, open(out_file, 'w') as o:
    for read in f:
        name = read.query_name
        cb = read.get_tag('XC')
        if cell.get(cb, 1) == 1:
            continue
        umi = read.get_tag('XM')
        all_tags = read.get_tags()
        for tags in all_tags:
            if 'GE' in tags:
                gene = tags[1]
                o.write(name + '\t' + cb + umi + '\t' + gene + '\n')
                break