#!/usr/bin/env python
# -*- coding: utf-8 -*-
# convert_primateai2tsi.py
# made by Daniel Minseok Kwon
# 2020-01-28 09:42:38
#########################
import sys
import os
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path = "/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path = "/ms1/bin/python_lib"
else:
    sys_path = "/home/mk446/bin/python_lib"
sys.path.append(sys_path)


def convert_primateai2tsi(orifile, out):
    
    header = ''
    f = {}
    for line in file_util.gzopen(orifile):
        line = line.decode('UTF-8')
        if line[0] != '#' and len(line.strip()) > 3:
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            if arr[0] == "chr":
                fields = '|'.join(arr[4:])
                arr[0] = "#CHROM"
                arr[1] = "POS"
                arr[2] = "ID"
                arr[3] = "REF"
                arr[4] = "ALT"
                arr[5] = fields
                header = '\t'.join(arr[:6]) + '\n'
            else:
                chrom = arr[0].replace('chr','')
                if not chrom in f.keys():
                    f[chrom] = open(out.replace('#CHROM#', chrom), 'w')
                    f[chrom].write(header)
                fields = '|'.join(arr[4:])
                cont = [chrom]
                cont.append(arr[1])
                cont.append('')
                cont.append(arr[2])
                cont.append(arr[3])
                cont.append(fields)
                f[chrom].write('\t'.join(cont) + '\n')
    for chrom in f.keys():
        f[chrom].close()

if __name__ == "__main__":
    import proc_util
    import file_util
    path = "/home/mk446/mutanno/DATASOURCE/PATHOGENICITY/PrimateAI/hg38/"
    orifile = path + "PrimateAI_scores_v0.2_hg38.tsv.gz"
    out = path + "PrimateAI_scores.v0.2.hg38.#CHROM#.tsi"
    convert_primateai2tsi(orifile, out)
