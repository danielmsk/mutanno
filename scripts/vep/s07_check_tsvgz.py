#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### s07_check_tsvgz.py
#### made by Daniel Minseok Kwon
#### 2020-01-20 16:23:31
#########################
import sys
import os
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path="/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path="/ms1/bin/python_lib"
else:
    sys_path="/home/mk446/bin/python_lib"
sys.path.append(sys_path)

import file_util
import proc_util


def s07_check_tsvgz():
    for fname in file_util.walk(path, 'tsv.vcf.gz'):
        # print(fname)
        if not file_util.is_exist(fname + '.tbi'):
            print(fname)

def check_chromvcf(chromvcf):
    i = 0
    for line in file_util.gzopen(chromvcf):
        line = line.decode('UTF-8')
        i += 1
        # if i >= 468687200 and i % 100000:
        #     print(i, line)
        # if '' in line:


# (base) mk446@compute-e-16-237:~/mutanno/PRECALVEP$ tabix -f -p vcf vep.98.hg38.1.tsv.gz
# [ti_index_core] the file out of order at line 468687254

if __name__ == "__main__":
    path = "/home/mk446/mutanno/PRECALVEP/chr1"
    chromvcf = "/home/mk446/mutanno/PRECALVEP/vep.98.hg38.1.tsv.gz"
    # s07_check_tsvgz()
    check_chromvcf(chromvcf)
