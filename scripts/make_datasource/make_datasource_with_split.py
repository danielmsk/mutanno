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


def merge():
    # f = open(out, 'w')
    logfile = out + '.log'
    chunklist = seq_util.get_split_region(100000, 'b38d')
    # print('total chunk size:', len(chunklist))
    flag = True
    for c1 in chunklist:
        tsvgz = path + str_util.zero_format(c1[3], 5) + ".tsv.gz"
        file_util.fileSave(logfile, tsvgz + '\n', 'a')
        for line in file_util.gzopen(tsvgz):
            line = line.decode('UTF-8')
            if line[0] != '#' or flag:
                # f.write(line)
                print(line, end='')
                flag = False
    # f.close()


def split_run():
    chunklist = seq_util.get_split_region(100000, 'b38d')
    print('total chunk size:', len(chunklist))
    for c1 in chunklist:
        # cmd = "python /home/mk446/mutanno/SRC/scripts/make_datasource/make_datasource_with_split.py "
        cmd = "mutanno makedata -ds /home/mk446/mutanno/SRC/tests/datastructure_microannot_v1.0.json"
        cmd += " -out " + path + str_util.zero_format(c1[3], 5) + ".tsv"
        cmd += " -region " + str(c1[0]) + ":" + str(c1[1]) + "-" + str(c1[2])
        print(cmd)


if __name__ == "__main__":
    import seq_util
    import str_util
    import file_util
    path = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/microanot_datasource_v1.0_tmp/"
    out = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/microanot_datasource_v1.0.tsv"
    if len(sys.argv) == 1:
        split_run()
    else:
        merge()
        # python make_datasource_with_split.py merge | bgzip -c > /home/mk446/mutanno/DATASOURCE/MICROANNOT/microanot_datasource_v1.0.tsv.gz
