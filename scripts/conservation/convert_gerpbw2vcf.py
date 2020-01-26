#!/usr/bin/env python
# -*- coding: utf-8 -*-
# convert_gerpbw2vcf.py
# made by Daniel Minseok Kwon
# 2020-01-25 22:56:45
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


def convert_gerpbw2vcf(bwfile):
    seq_util.load_refseq_info('b38d')
    for chrom in seq_util.MAIN_CHROM_LIST:
        if chrom == "MT":
            chrom = "M"
        cmd = "bigWigToWig " + bwfile + ' ' + out.replace('#CHROM#', chrom)
        cmd += ' -chrom=' + chrom
        # cmd += ' -start=1'
        # cmd += ' -end=' + str(seq_util.CHROM_LEN['b38d'][chrom] + 1)
        print(cmd)


if __name__ == "__main__":
    import seq_util
    bwfile = "/home/mk446/mutanno/DATASOURCE/CONSERVATION/GERP/hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
    out = bwfile + '.#CHROM#.wig'
    convert_gerpbw2vcf(bwfile)
