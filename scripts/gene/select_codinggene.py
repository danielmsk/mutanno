#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ./scripts/gene/select_codinggene.py
# made by Daniel Minseok Kwon
# 2020-03-31 04:09:37
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


def select_codinggene(infile, outfile):
    f = open(outfile, 'w')
    h = {}
    for line in file_util.gzopen(infile):
        line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        flag = False
        if line[0] == "#":
            flag = True
            for k in range(len(arr)):
                h[arr[k]] = k
        else:
            # cidx = h['HGNC_locus_group']
            # if arr[cidx] == "protein-coding gene":
            #     flag = True

            cidx = h['ENSEMBLgene_gene_biotype']
            if arr[cidx] == "protein_coding":
                flag = True

        if flag:
            f.write(line)

    f.close()   

if __name__ == "__main__":
    import proc_util
    import file_util
    infile = sys.argv[1]
    outfile = sys.argv[2]
    select_codinggene(infile, outfile)
