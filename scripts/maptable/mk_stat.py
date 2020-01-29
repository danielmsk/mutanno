#!/usr/bin/env python
# -*- coding: utf-8 -*-
# mk_stat.py
# made by Daniel Minseok Kwon
# 2020-01-29 09:41:26
#########################
import sys
import os
import json
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path = "/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path = "/ms1/bin/python_lib"
else:
    sys_path = "/home/mk446/bin/python_lib"
sys.path.append(sys_path)


def mk_stat(jsonfile):
    ds = ""
    with open(jsonfile) as jfp:
        ds = json.load(jfp)

    for s1 in ds['source']:
        print(os.path.join(ds['datafile_path'], s1['datafile']))

if __name__ == "__main__":
    import proc_util
    import file_util
    jsonfile = "/home/mk446/mutanno/SRC/tests/datastructure_v0.3.0_mvp.json"
    mk_stat(jsonfile)
