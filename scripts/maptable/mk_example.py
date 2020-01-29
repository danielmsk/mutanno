#!/usr/bin/env python
# -*- coding: utf-8 -*-
# mk_example.py
# made by Daniel Minseok Kwon
# 2020-01-29 15:06:07
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

def load_clinvar():
    m = []
    for line in file_util.gzopen(clinvarfile):
        line = line.decode('UTF-8')
        if line[0] != "#":
            arr = line.split('\t')
            if arr[0] != '1':
                break
            if arr[11] == "Pathogenic" and len(arr[2]) == 1 and len(arr[3]) == 1:
                m.append(arr[:4])
    return m

def mk_example():
    m = load_clinvar()
    f = open(out, 'w')
    k = 0
    for line in file_util.gzopen(base_vcf):
        line = line.decode('UTF-8')
        if line[0] == "#":
            f.write(line)
        else:
            arr = line.split('\t')
            if arr[6] == "PASS" and len(arr[3]) == 1 and len(arr[4]) == 1:
                arr[1] = m[k][1]
                arr[2] = '.'
                arr[3] = m[k][2]
                arr[4] = m[k][3]
                print (arr)
                f.write('\t'.join(arr))
                k += 1
                if k >= len(m):
                    break
    f.close()
    print('Saved',out)

if __name__ == "__main__":
    import proc_util
    import file_util
    base_vcf = "/home/mk446/bio/DATA/AshkenazimTrio/hg38/TRIO.x60/Trio.hs38d1.60x.1.GVCF.mnv.vqsr.vcf.gz"
    out = "/home/mk446/mutanno/DATASOURCE/TEST/test_trio_v0.3.vcf"
    clinvarfile = "/home/mk446/mutanno/DATASOURCE/VARIANTDB/CLINVAR/hg38/clinvar.20200106.hg38.tsv.gz"
    mk_example()
