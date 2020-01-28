#!/usr/bin/env python
# -*- coding: utf-8 -*-
# split_20way_by_chrom.py
# made by Daniel Minseok Kwon
# 2020-01-27 14:42:24
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


def split_20way_by_chrom(wigfile, outwig):
    f = ''
    prev_chrom = ''
    for line in file_util.gzopen(wigfile):
        line = line.decode('UTF-8')
        if line[:len('fixedStep')] == "fixedStep":
            arr = line.split(' ')
            arr[-1] = arr[-1].strip()
            chrom = arr[1].replace('chrom=chr', '')
            if chrom != prev_chrom:
                try:
                    f.close()
                except AttributeError:
                    pass
                out = outwig.replace('#CHROM#', chrom)
                print(out)
                f = open(out, 'w')
            print(line.strip(), out)
            prev_chrom = chrom
        f.write(line)
    f.close()


if __name__ == "__main__":
    import proc_util
    import file_util
    wigfile = "/home/mk446/mutanno/DATASOURCE/CONSERVATION/PHASTCONS/PHASTCONS20WAY/hg38/hg38.phastCons20way.wigFix.gz"
    outwig = "/home/mk446/mutanno/DATASOURCE/CONSERVATION/PHASTCONS/PHASTCONS20WAY/hg38/chr#CHROM#.phastCons20way.wigFix"
    # split_20way_by_chrom(wigfile, outwig)
    wigfile = "/home/mk446/mutanno/DATASOURCE/CONSERVATION/PHYLOP/PHYLOP20WAY/hg38/hg38.phyloP20way.wigFix.gz"
    outwig = "/home/mk446/mutanno/DATASOURCE/CONSERVATION/PHYLOP/PHYLOP20WAY/hg38/chr#CHROM#.phyloP20way.wigFix"
    split_20way_by_chrom(wigfile, outwig)
