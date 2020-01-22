#!/usr/bin/env python
# -*- coding: utf-8 -*-
# rm_redundant_variant.py
# made by Min-Seok Kwon
# 2020-01-14 11:16:03
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


def rm_redundant_variant(vcf, out):
    prev_line = ''
    f = open(out, 'w')
    for line in open(vcf):
        if line != prev_line:
            f.write(prev_line)
        prev_line = line
    f.write(prev_line)
    f.close()


if __name__ == "__main__":
    vcf = '/home/mk446/mutanno/DATASOURCE/KNOWN_INDEL/hg38/known_indel.sorted.vcf'
    out = '/home/mk446/mutanno/DATASOURCE/KNOWN_INDEL/hg38/known_indel.sorted.uniq.vcf'
    rm_redundant_variant(vcf, out)
