#!/usr/bin/env python
# -*- coding: utf-8 -*-
# merge_datasource.py
# made by Daniel Minseok Kwon
# 2020-01-29 14:24:50
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

def zero_format(i, size):
    # 00001000 format
    s1 = "{:0>" + str(size) + "d}"
    return s1.format(i)


def merge_datasource():
    pos = 1
    k = 0
    flag = True
    while flag:
        k += 1
        epos = pos + bsize - 1
        if epos >= seq_util.CHROM_LEN['b38d']['1']:
            epos = seq_util.CHROM_LEN['b38d']['1']
            flag = False
        
        sfile = path + zero_format(k, 5) + '.tsi'
        if k == 1:
            cmd = "cat " + sfile + " > " + out + ";" 
        else:
            cmd = "cat " + sfile + " | grep -v '^#' >> " + out + ";"
        cmd += "echo '" + sfile + "';"
        print(cmd)
        pos = epos - 1
    
    print ('sleep 120;')
    print ('tabixgz ' + out)
    

if __name__ == "__main__":
    import seq_util
    import file_util
    seq_util.load_refseq_info('b38d')
    bsize = 100000
    out = "/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/mvp_datasource_v0.3.chr1.tsi"
    path = "/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/mvp_datasource_v0.3_test_tmp/"
    merge_datasource()
