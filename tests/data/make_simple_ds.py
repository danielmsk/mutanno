#!/usr/bin/env python
# -*- coding: utf-8 -*-
# make_simple_ds.py
# made by Daniel Minseok Kwon
# 2020-04-30 11:20:00
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
import json


def make_simple_ds(dsfile):
    out = dsfile.replace('ds.json', '.json') + '.simple'

    ds = file_util.jsonOpen(dsfile)

    dskeylist = list(ds.keys())
    for k1 in dskeylist:
        if k1 in ['datafile_path','sourcefile','use_sourcename_in_fieldname']:
            del(ds[k1])
        # print(ds[k1])

    for i in range(len(ds['source'])):
        dskeylist = list(ds['source'][i].keys())
        for k1 in dskeylist:
            if k1 in ['datafile','format','ref_column_index','alt_column_index','fieldselection']:
                del(ds['source'][i][k1])

            del_fieldname_list = []
            fieldlen = len(ds['source'][i]['fields'])
            for j in range(fieldlen):
                if 'is_available' in ds['source'][i]['fields'][j].keys():
                    if not ds['source'][i]['fields'][j]['is_available']:
                        fn = ds['source'][i]['fields'][j]['name']
                        del_fieldname_list.append(fn)

            print(del_fieldname_list)
            for fn in del_fieldname_list:
                fieldlen = len(ds['source'][i]['fields'])
                for j in range(fieldlen):
                    if fn == ds['source'][i]['fields'][j]['name']:
                        del(ds['source'][i]['fields'][j])


    file_util.jsonSave(out, ds, 2)
    print('Saved', out)

if __name__ == "__main__":
    import proc_util
    import file_util
    # dsfile = "datastructure_microannot_v0.4.1ds.json"
    dsfile = sys.argv[1]
    make_simple_ds(dsfile)
