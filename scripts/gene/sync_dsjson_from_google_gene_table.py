#!/usr/bin/env python
# -*- coding: utf-8 -*-
# sync_field_avail.py
# made by Daniel Minseok Kwon
# 2020-04-03 16:15:24
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

whitefields = ['ensgid']

def check_two_fields(json_dict, json_key,  maptab_dict, maptab_key, fn):
    if json_key in json_dict.keys():
        if json_dict[json_key] and maptab_dict[maptab_key].upper() == 'Y':
            pass
        elif not json_dict[json_key] and maptab_dict[maptab_key].upper() == 'N':
            pass
        else:
            print('\t', json_key, json_dict[json_key], '=>', maptab_dict[maptab_key].upper())
            json_dict[json_key] = maptab_dict[maptab_key].upper() == 'Y'
    else:
        print('\t', json_key, fn, 'NA =>', maptab_dict[maptab_key].upper())
        json_dict[json_key] = maptab_dict[maptab_key].upper() == 'Y'

    return json_dict

def check_two_fields2(json_dict, json_key,  maptab_dict, maptab_key, fn):
    if json_key in json_dict.keys():
        if json_dict[json_key] != maptab_dict[maptab_key]:
            print('\t', json_key, json_dict[json_key], '=>', maptab_dict[maptab_key])
            json_dict[json_key] = maptab_dict[maptab_key]
    else:
        print('\t', json_key, fn, 'NA =>', maptab_dict[maptab_key])
        json_dict[json_key] = maptab_dict[maptab_key]

    return json_dict

def get_fieldname(field):
    if 'name2' in field.keys():
        fn = field['name2']
    else:
        fn = field['name']
    return fn


def sync_field_avail(json_file, maptab_file, json_out):


    maptab = {}
    h = {}
    for line in open(maptab_file):
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        if arr[0] == "no":
            # print(arr)
            for k in range(len(arr)):
                h[arr[k]] = k
        else:
            # print(arr[h['Datasource']], arr[h['Field Name']], arr[h['do import']], arr[h['do_import_2']])
            ds = arr[h['Datasource']]
            try:
                maptab[ds]
            except KeyError:
                maptab[ds] = {}

            fn = arr[h['Field Name']]

            d = {}
            for k1 in h.keys():
                d[k1] = arr[h[k1]]
            maptab[ds][fn] = d


    with open(json_file) as f:
        data = json.load(f)

    # print(data['gene_source'].keys())

    for gs in data['gene_source']:
        # print(gs.keys())
        print('DS:',gs['name'])

        ds = gs['name']
        maptab[ds]
        for field in gs['fields']:
            fn = get_fieldname(field)
            if not fn in whitefields:
                if fn in maptab[ds].keys():
                    maptabdata = maptab[ds][fn]
                    field = check_two_fields(field,'do_import', maptabdata, 'do_import_2', fn)
                    field = check_two_fields(field,'is_list', maptabdata, 'is_list', fn)
                    field = check_two_fields2(field,'type', maptabdata, "field_type (string, integer, boolean, number)", fn)
                    
                    if 'is_available' in field.keys() and not field['is_available']:
                        print('\t', 'is_available',fn, field['is_available'], '=>', True)
                        field['is_available'] = True
                else:
                    if field['is_available']:
                        print('\t', 'is_available', fn, field['is_available'], '=>', False)
                        field['is_available'] = False

        # break


    with open(json_out, 'w') as f2:
        json.dump(data, f2, indent=4)

    print('Saved', json_out)
    

if __name__ == "__main__":
    import proc_util
    import file_util
    json_file = "/home/mk446/mutanno/SRC/tests/data/datastructure_v0.4.0ds_mvp.json"
    ### from https://docs.google.com/spreadsheets/d/1HYTkkhTcW_rsaYBcjjL-I4tgLlO9o0Yuz1rCvifwMRs/edit#gid=853902134
    maptab_file = "/home/mk446/mutanno/SRC/tests/data/datastructure_v0.4.0ds_mvp.vcfmaptab.txt"

    json_out = "/home/mk446/mutanno/SRC/tests/data/datastructure_v0.4.0ds_mvp_out.json"
    
    sync_field_avail(json_file, maptab_file, json_out)
