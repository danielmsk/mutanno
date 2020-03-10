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
import time


def merge(chrom):
    out2 = out.replace('#CHROM#', chrom)
    f = open(out2, 'w')
    print('Saving..',out2)
    chunklist = seq_util.get_split_region(binsize, 'b38d')
    # print('total chunk size:', len(chunklist))
    flag_header = True
    for c1 in chunklist:
        # if chrom != c1[0]:
        #     f.close()
        #     chrom = c1[0]
        #     out2 = out.replace('#CHROM#', chrom)
        #     f = open(out2, 'w')
        #     print('Saving..',out2)
        #     flag_header = True
        if chrom == c1[0]:
            tsvgz = tmp_path + str_util.zero_format(c1[3], 5) + ".tsi.gz"
            print (tsvgz)
            for line in file_util.gzopen(tsvgz):
                line = line.decode('UTF-8')
                if line[0] == '#' and flag_header:
                    f.write(line)
                    flag_header = False
                if line[0] != '#':
                    f.write(line)
    f.close()
    time.sleep(5)
    proc_util.run_cmd('tabixgz ' + out2)

def run_bychrom():
    for chrom in seq_util.MAIN_CHROM_LIST:
        print ("python " + os.path.abspath("./s03_merge_datasource.py") + " " + chrom)

if __name__ == "__main__":
    import seq_util
    import str_util
    import file_util
    import proc_util
    path = "../../../DATASOURCE/MUTANOANNOT/"
    # out = path + "microanot_datasource_v1.0.tsi"
    out = path + "mvp_datasource_v0.3.2.chr#CHROM#.tsi"
    tmp_path = path + "datasource_v0.3_tmp/"
    binsize = 1000000
    
    if len(sys.argv) == 1:
        run_bychrom()
    else:
        merge(sys.argv[1])
    # python make_datasource_with_split.py merge | bgzip -c > /home/mk446/mutanno/DATASOURCE/MICROANNOT/microanot_datasource_v1.0.tsv.gz
