#!/usr/bin/env python
# -*- coding: utf-8 -*-
# make_datasource_with_split.py
# made by Min-Seok Kwon
# 2020-01-17 16:40:44
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
import tabix


def check_file(chunk_file, region, vep):
    flag = True
    varmap = {}
    for line in file_util.gzopen(chunk_file):
        line = line.decode('UTF-8')
        if line[0] != "#":
            r1 = line.split('\t')
            k1 = r1[0] + '\t' + r1[1] + '\t' + r1[3]+ '\t' + r1[4]
            varmap[k1] = 0

    errmsg = ""
    tp = tabix.open(vep)
    for r1 in tp.querys(region):
        k1 = r1[0] + '\t' + r1[1] + '\t' + r1[3]+ '\t' + r1[4]
        try:
            tmp = varmap[k1]
        except KeyError:
            flag = False
            errmsg += k1 + '\n'
        

    if flag:
        file_util.fileSave(chunk_file + '.checked', '', 'w')
        print ('checked')
    else:
        file_util.fileSave(chunk_file + '.error', errmsg, 'w')
        print ('error')


def run_all(vep):
    chunklist = seq_util.get_split_region(binsize, 'b38d')
    for c1 in chunklist:
        chunk_file = os.path.abspath(tmp_path + str_util.zero_format(c1[3], 5) + ".tsi.gz")
        region = str(c1[0]) + ":" + str(c1[1]) + "-" + str(c1[2])
        # total_line = vep_tp.querys(region)
        cmd = "python " + os.path.abspath('s02_check_split_datasourcefiles.py') + " " + chunk_file 
        cmd += " " + region + " " + vep.replace('#CHROM#', c1[0])
        print(cmd)
    # print('#total chunk size:', len(chunklist))

if __name__ == "__main__":
    import seq_util
    import str_util
    import file_util
    path = "../../../DATASOURCE/MUTANOANNOT/"
    out = path + "microanot_datasource_v1.0.tsi"
    tmp_path = path + "datasource_v0.3_tmp/"
    dsfile = "../../tests/datastructure_v0.3.1_mvp.json"
    binsize = 1000000
    if len(sys.argv) == 1:
        vep = "/home/mk446/mutanno/DATASOURCE/ANNOT/VEP/hg38/v99/vep.99.hg38.#CHROM#.sorted.rmcsq.tsi.gz"
        run_all(vep)
    else:
        chunk_file = sys.argv[1]
        region = sys.argv[2]
        vep = sys.argv[3]
        check_file(chunk_file, region, vep)
    
