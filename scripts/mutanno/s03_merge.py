#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s03_merge.py
# made by Daniel Minseok Kwon
# 2020-02-11 10:31:49
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


def s03_merge(tchrom):
    regionlist = seq_util.get_split_region()
    
    flagheader = True
    prev_pos = 0
    out2 = path + out.replace('#CHROM#', tchrom)
    f = open(out2, 'w')
    for r1 in regionlist:
        chrom = r1[0]
        if chrom == tchrom:
            tsifile = path + "tmp/mc_" + str(r1[3]) + ".tsv.tsi"
            for line in open(tsifile):
                if line[0]=='#':
                    if flagheader:
                        f.write(line)
                        flagheader = False
                else:
                    arr = line.split('\t')
                    pos = int(arr[1])
                    if prev_pos <= pos:
                        f.write(line)
                        prev_pos = pos
                    else:
                        print('ERROR position', prev_pos, arr)
    f.close()
    time.sleep(30)
    proc_util.run_cmd('tabixgz ' + out2)

def run():
    regionlist = seq_util.get_split_region()
    prev_chrom = ''
    for r1 in regionlist:
        chrom = r1[0]
        if chrom != prev_chrom:
            cmd = "python /home/mk446/mutanno/SRC/scripts/microannotation/s03_merge.py " + chrom
            print(cmd)
        prev_chrom = chrom


if __name__ == "__main__":
    import proc_util
    import file_util
    import seq_util
    path = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/"
    out = "microannot_datasource.#CHROM#.v0.2_200211.tsi"
    if len(sys.argv) == 1:
        run()
    else:
        s03_merge(sys.argv[1])
