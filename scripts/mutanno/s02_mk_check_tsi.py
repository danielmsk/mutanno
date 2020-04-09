#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s02_mk_check_tsi.py
# made by Daniel Minseok Kwon
# 2020-02-10 09:04:55
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


def s02_mk_check_tsi(tsifile, spos, epos):
    flag = True
    tpos = spos
    pos_flag = False
    for line in open(tsifile):
        if line[0] != '#':
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            if "VEP=" not in arr[-1]:
                flag = False
                print(arr)
            # pos = int(arr[1])
            # if pos == tpos + 1 and pos_flag:
            #     pos_flag = True
            #     tpos = pos + 1
            # elif pos == tpos:
            #     pos_flag = True
            # else:
            #     pos_flag = False
            #     print("ERROR ", arr[0], arr[1], tpos)
            #     break
            # if not pos_flag:
            #     flag = False

    if flag:
        donefile = tsifile + ".done"
        print(donefile)
        file_util.fileSave(donefile, '', 'w')
    

def run():
    regionlist = seq_util.get_split_region()
    for r1 in regionlist:
        tsifile = path + "mc_" + str(r1[3]) + ".tsv.tsi"
        cmd = "python /home/mk446/mutanno/SRC/scripts/mutanno/s02_mk_check_tsi.py " + tsifile
        cmd += " " + str(r1[1])
        cmd += " " + str(r1[2])
        print(cmd)


if __name__ == "__main__":
    import proc_util
    import file_util
    import seq_util
    path = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/tmp/"
    if len(sys.argv) == 1:
        run()
    else:
        tsifile = sys.argv[1]
        spos = int(sys.argv[2])
        epos = int(sys.argv[3])
        s02_mk_check_tsi(tsifile, spos, epos)
