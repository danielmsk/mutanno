#!/usr/bin/env python
# -*- coding: utf-8 -*-
# convert_transmem.py
# made by Daniel Minseok Kwon
# 2020-01-27 16:23:07
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


def convert_transmem():
    out = bedfile.replace('.bed.gz', '.mod.bed')
    f = open(out, 'w')
    f.write('#CHROM\tSPOS\tEPOS\tUNIPROTKB_AC\tSCORE\tSTRAND\tTHICK_START\tTHICK_END\tANNOTATION_COLOR\tNO_BLOCK\tBLOCK_SIZE\tBLOCK_START\tANNOT_ID\tANNOT_DESC\n')
    for line in file_util.gzopen(bedfile):
        line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[0] = arr[0].replace('chr', '')
        if arr[0] == "MT":
            arr[0] = "M"
        f.write('\t'.join(arr))
    f.close()
    proc_util.run_cmd('tabixgzbed ' + out)


if __name__ == "__main__":
    import proc_util
    import file_util
    bedfile = "/home/mk446/mutanno/DATASOURCE/UNIPLOT/UP000005640_9606_transmem.bed.gz"
    convert_transmem()
