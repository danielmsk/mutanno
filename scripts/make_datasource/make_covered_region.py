#!/usr/bin/env python
# -*- coding: utf-8 -*-
# make_covered_region.py
# made by Daniel Minseok Kwon
# 2020-05-05 04:39:33
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


def make_covered_region(vcf, outbed):
    f = open(outbed, 'w')

    spos = -9
    ipos = -9
    for line in file_util.gzopen(vcf):
        line = file_util.decodeb(line)
        if line[0] != '#':
            arr = line.split('\t')
            chrom = arr[0]
            pos = int(arr[1])

            if ipos == pos:
                pass
            elif ipos + 1 == pos:
                ipos += 1
            elif ipos + 1 < pos and ipos > 0:
                epos = ipos
                cont = chrom + '\t' + str(spos) + '\t' + str(epos) + '\n'
                f.write(cont)
                spos = pos -1 
                ipos = pos

            if spos < 0:
                spos = pos - 1
                ipos = pos

    if spos + 1 < ipos:
        epos = ipos
        cont = chrom + '\t' + str(spos) + '\t' + str(epos) + '\n'
        f.write(cont)

    f.close()
    

if __name__ == "__main__":
    import proc_util
    import file_util
    vcf = sys.argv[1]
    outbed = vcf + '.covered.bed'
    make_covered_region(vcf, outbed)
