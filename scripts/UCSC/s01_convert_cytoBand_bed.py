#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s01_convert_cytoBand_bed2vcf.py
# made by Daniel Minseok Kwon
# 2020-01-27 11:10:27
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


def s01_convert_cytoBand_bed2vcf(bed, out):
    f = open(out, 'w')
    header = "#CHROM\tSPOS\tEPOS\tcytoband"
    f.write(header + '\n')
    for line in file_util.gzopen(bed):
        line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        chrom = arr[0].replace('chr', '')
        cont = [chrom]
        cont.append(arr[1])
        cont.append(arr[2])
        cont.append(chrom + arr[3])
        f.write('\t'.join(cont) + '\n')
    f.close()


if __name__ == "__main__":
    import file_util
    bed = "/home/mk446/mutanno/DATASOURCE/CYTOBAND/hg38/cytoBand.txt.gz"
    out = "/home/mk446/mutanno/DATASOURCE/CYTOBAND/hg38/cytoBand_hg38.bed"
    s01_convert_cytoBand_bed2vcf(bed, out)
