#!/usr/bin/env python
# -*- coding: utf-8 -*-
# reformat.py
# made by Daniel Minseok Kwon
# 2020-01-28 04:36:09
#########################
import sys
import os
import time
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path = "/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path = "/ms1/bin/python_lib"
else:
    sys_path = "/home/mk446/bin/python_lib"
sys.path.append(sys_path)


def reformat(tsi):
    out = path + "tmp2/" + tsi.split('/')[-1].replace('.tsi.gz', '') + '.bed'
    f = open(out, 'w')
    for line in file_util.gzopen(tsi):
        if tsi.endswith('.gz'):
            line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        if line[0] == '#':
            spos = "SPOS"
            epos = "EPOS"
        else:
            spos = str(int(arr[1]) - 1)
            epos = arr[1]
        cont = [arr[0], spos, epos, arr[3].replace('|','\t')]
        f.write('\t'.join(cont) + '\n')
    f.close()
    print("Saved", out)

    time.sleep(120)
    
    proc_util.run_cmd('tabixgzbed ' + out)


# /home/mk446/mutanno/SRC/scripts/conservation/reformat.py

if __name__ == "__main__":
    import file_util
    import proc_util
    tsi = sys.argv[1]
    path = "/home/mk446/mutanno/DATASOURCE/CONSERVATION/"
    # tsi = path + "conservation_scores.hg38.chrM.tsi"
    reformat(tsi)
