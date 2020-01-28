#!/usr/bin/env python
# -*- coding: utf-8 -*-
# split_variants_keep_only_id.py
# made by Daniel Minseok Kwon
# 2020-01-27 17:17:25
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


def split_variants_keep_only_id():
    f = open(out, 'w')
    for line in file_util.gzopen(dbsnpfile):
        line = line.decode('UTF-8')
        if line[0] == '#':
            if line[:len('#CHROM')] == "#CHROM":
                arr = line.split("\t")
                f.write('\t'.join(arr[:5]) + '\n')
        else:
            arr = line.split('\t')
            for alt in arr[4].split(','):
                f.write('\t'.join([arr[0], arr[1], arr[2], arr[3], alt]) + '\n')
    f.close()


if __name__ == "__main__":
    import proc_util
    import file_util
    dbsnpfile = "/home/mk446/mutanno/DATASOURCE/dbSNP/hg38/GRCh38_latest_dbSNP_all.convchrom.vcf.gz"
    out = dbsnpfile.replace('.vcf.gz', '') + '.onlyid.vcf'
    split_variants_keep_only_id()
