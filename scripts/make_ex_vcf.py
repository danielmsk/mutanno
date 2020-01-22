#!/usr/bin/env python
# -*- coding: utf-8 -*-
# make_ex_vcf.py
# made by Min-Seok Kwon
# 2020-01-14 15:36:32
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


def make_ex_vcf(clinvar_vcf, real_vcf, out, no):
    f = open(out, 'w')
    arrvar = []
    i = 0
    for line in file_util.gzopen(real_vcf):
        line = line.decode('UTF-8')
        if line[0] == '#':
            f.write(line)
        else:
            arr = line.split('\t')
            arrvar.append(arr)
            i += 1
            if i > no:
                break

    i = 0
    for line in file_util.gzopen(clinvar_vcf):
        line = line.decode('UTF-8')
        if line[0] != '#':
            arr = line.split('\t')
            arr2 = arrvar[i]
            arr2[0] = arr[0]
            arr2[1] = arr[1]
            arr2[2] = arr[2]
            arr2[3] = arr[3]
            arr2[4] = arr[4]
            f.write('\t'.join(arr2))
            i += 1
            if i >= len(arrvar):
                break
    f.close()


if __name__ == "__main__":
    import proc_util
    import file_util
    clinvar_vcf = "/home/mk446/mutanno/DATASOURCE/VARIANTDB/CLINVAR/hg38/clinvar.20200106.hg38.tsv.gz"
    real_vcf = "/home/mk446/bio/DATA/AshkenazimTrio/hg38/TRIO.x60/Trio.hs38d1.60x.1.GVCF.mnvindel.vqsr.vcf.gz"
    out = "/home/mk446/mutanno/DATASOURCE/TEST/trio_clinvar_variants_100.vcf"
    make_ex_vcf(clinvar_vcf, real_vcf, out, 100)
