#!/usr/bin/env python
# -*- coding: utf-8 -*-
# tmp_add_na_microannot.py
# made by Daniel Minseok Kwon
# 2020-03-05 11:40:02
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
import time

def tmp_add_na_microannot(infile, outfile):
    chrom = infile.split('/')[-1].split('.')[1]
    out = outfile.replace('#CHROM#',chrom)
    f = open(out, 'w')
    i = 0
    for line in file_util.gzopen(infile):
        line = line.decode('UTF-8')
        if line[0] == "#":
            f.write(line)
        else:
            i += 1
            if "gnomADgenome=" in line:
                f.write(line)
            else:
                arr = line.split('\t')
                arr[-1] = arr[-1].strip()
                if arr[-1] != "":
                    arr[-1] += ";"    
                arr[-1] += "gnomADgenome=NA"
                # print(arr[-1])
                f.write('\t'.join(arr) + '\n')
                # if i > 10:
                #     break
    f.close()
    time.sleep(10)
    proc_util.run_cmd('tabixgz ' + out)


def run(path):
    for fname in file_util.walk(path, '.tsi.gz'):
        cmd = "python tmp_add_na_microannot.py " + fname
        print(cmd)

if __name__ == "__main__":
    import proc_util
    import file_util
    path = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/v0.2_200211/"
    outpath = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/v0.3_200305/"
    outfile = outpath + "microannot_datasource.#CHROM#.v0.3_200305.tsi"
    if len(sys.argv) == 1:
        run(path)
    else:
        infile = sys.argv[1]
        tmp_add_na_microannot(infile, outfile)
