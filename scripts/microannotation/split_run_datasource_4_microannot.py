#!/usr/bin/env python
# -*- coding: utf-8 -*-
# split_run_datasource_4_microannot.py
# made by Min-Seok Kwon
# 2020-01-14 12:14:40
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


def split_run_datasource_4_microannot():
    # seq_util.load_refseq_info('b38d')
    regionlist = seq_util.get_split_region()
    # print(regionlist[:20])
    # print(len(regionlist))
    for r1 in regionlist:
        cmd = 'python /home/mk446/mutanno/SRC/mutanno.py makedata '
        cmd += "-ds /home/mk446/mutanno/SRC/tests/datastructure_microannot_v1.0.json "
        cmd += "-out /home/mk446/mutanno/DATASOURCE/MICROANNOT/tmp/" + str(r1[3]) + ".tsv "
        cmd += "-region " + r1[0] + ":" + str(r1[1]) + "-" + str(r1[2]) + ";"
        # print(cmd)
        sh = path + 'aaaaa_' + str(r1[3]) + ".tsv.sh"
        file_util.fileSave(sh, cmd + '\n', 'w')
    proc_util.run_cmd('chmod 755 ' + path + '*.sh')


if __name__ == "__main__":
    import proc_util
    import file_util
    import seq_util
    path = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/tmp/"
    split_run_datasource_4_microannot()
