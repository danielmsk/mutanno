#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### s09_check_chrom1.py
#### made by Min-Seok Kwon
#### 2020-01-21 09:55:02
#########################
import sys
import os
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path="/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path="/ms1/bin/python_lib"
else:
    sys_path="/home/mk446/bin/python_lib"
sys.path.append(sys_path)
import file_util
import proc_util

def s09_check_chrom1(chrom):
    i = 0
    for k in range(400):
        vcfmap = {}
        for tsv in file_util.walk(path + "chr" + chrom + "/" + str(k) + "/", '_microannot.tsv'):
            k1 = int(tsv.split('/')[-1].split('_')[1])
            vcfmap[tsv] = k1
        (ks, vs) = struct_util.sortdict(vcfmap)

        prev_pos = 0
        for tsv in ks:
            j = 0
            log = str(i + 1) + ': ' + tsv
            print(log)
            for line in file_util.gzopen(tsv):
                # line = line.decode('UTF-8')
                if (i == 0 and line[0] == '#'):
                    line = line.replace('\t.\t', '\tID\t')
                if (i == 0) or (i > 0 and line[0] != '#'):
                    # print(line, end='')
                    if line[0] != '#':
                        arr = line.split('\t')
                        pos = int(arr[1])
                        if pos < prev_pos:
                            print('ERROR', i, arr, prev_pos)
                        prev_pos = pos
                j += 1
                if j % 100000 == 0:
                    print('\t' + line[:30])
            i += 1
        


if __name__ == "__main__":
    import struct_util
    path = "/home/mk446/mutanno/PRECALVEP/"
    s09_check_chrom1("1")
