#!/usr/bin/env python
# -*- coding: utf-8 -*-
# preproc_dbNSFP.py
# made by Daniel Minseok Kwon
# 2020-03-05 08:54:00
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


def preproc_dbNSFP(infile, outfile):
    f = open(outfile, 'w')
    i = 0
    for line in file_util.gzopen(infile):
        line = line.decode('UTF-8')
        arr = line.split('\t')
        if i == 0:
            f.write('#'+line)
        else:
            arr[-1] = arr[-1].strip()
            for k in range(len(arr)):
                arr[k] = arr[k].strip()
                if arr[k] == ".":
                    arr[k] = ""
            f.write('\t'.join(arr) + '\n')
        i += 1

    f.close()

    proc_util.run_cmd('tabixgzbed ' + outfile)
    

if __name__ == "__main__":
    import proc_util
    import file_util
    path = "/home/mk446/mutanno/DATASOURCE/VARIANTDB/dbNSFP/hg38/dbNSFP4.0c/"
    infile = path + "dbNSFP4.0_gene.complete.gz"
    outfile = path + "dbNSFP4.0_gene.complete.mod"
    preproc_dbNSFP(infile, outfile)
