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


def extract_indel():
    f = open(out, 'w')
    for line in file_util.gzopen(onlyidfile):
        line = line.decode('UTF-8')
        if line[0] == '#':
            if line[:len('#CHROM')] == "#CHROM":
                arr = line.split("\t")
                arr[-1] = arr[-1].strip()
                f.write('\t'.join(arr[:5]) + '\n')
        else:
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            for alt in arr[4].split(','):
                if len(arr[3]) > 1 or len(alt) > 1:
                    f.write('\t'.join([arr[0], arr[1], arr[2], arr[3], alt]) + '\n')
    f.close()
    print("Saved", out)


if __name__ == "__main__":
    import proc_util
    import file_util
    onlyidfile = "/home/mk446/bio/mutanno/DATASOURCE/dbSNP/hg38/GRCh38_latest_dbSNP_all.convchrom.onlyid.vcf.gz"
    out = onlyidfile.replace('.vcf.gz', '') + '.indel.vcf'
    extract_indel()
