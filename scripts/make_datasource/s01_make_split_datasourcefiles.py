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


def split_run():
    chunklist = seq_util.get_split_region(binsize, 'b38d')
    for c1 in chunklist:
        # cmd = "python /home/mk446/mutanno/SRC/scripts/make_datasource/make_datasource_with_split.py "
        tmpout = os.path.abspath(tmp_path + str_util.zero_format(c1[3], 5) + ".tsi")
        cmd = "mutanno makedata -ds " + os.path.abspath(dsfile)
        cmd += " -out " + tmpout
        cmd += " -region " + str(c1[0]) + ":" + str(c1[1]) + "-" + str(c1[2])
        cmd += ";"
        cmd += "sleep 3;"
        cmd += "tabixgz " + tmpout + ";"
        cmd += "sleep 3;"
        cmd += "rm " + tmpout + ";"
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
    split_run()
