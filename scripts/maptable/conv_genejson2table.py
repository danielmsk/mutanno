#!/usr/bin/env python
# -*- coding: utf-8 -*-
# conv_genejson2table.py
# made by Daniel Minseok Kwon
# 2020-02-04 13:40:44
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

def dv(d1, k1):
    v1 = ''
    if k1 in d1.keys():
        v1 = d1[k1]
    return v1

def load_filed_stat_from_datafile():
    stat = {}
    stat['maxlen'] = {}
    stat['cnt'] = {}
    header = []
    h = {}
    i = 0
    for line in file_util.gzopen(datafile):
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        if line[0] == "#":
            arr[0]= arr[0][1:]
            header = arr
            for j in range(len(arr)):
                h[arr[j]] = j
                
        else:
            i += 1
            for j in range(len(arr)):
                k1 = header[j]
                for xid in arr[js].split('|'):
                    if len(xid) > stat['maxlen'][k1]
                try:
                    if len(values[k]) > stat['maxlen'][k1]:
                        stat['maxlen'][k1] = len(values[k])
                except KeyError:
                    stat['maxlen'][k1] = len(values[k])

                sourcename = arr2[0]
                for section in arr2[1].split(','):
                    values = section.split('|')
                    for k in range(len(values)):
                        k1 = sourcename.lower() + '_' + h[sourcename][k].lower()
                        try:
                            if len(values[k]) > stat['maxlen'][k1]:
                                stat['maxlen'][k1] = len(values[k])
                        except KeyError:
                            stat['maxlen'][k1] = len(values[k])

                        try:
                            tmp = stat['cnt'][k1]
                        except KeyError:
                            stat['cnt'][k1] = {}

                        try:
                            stat['cnt'][k1][values[k]] += 1
                        except KeyError:
                            stat['cnt'][k1][values[k]] = 1
            # break
            if i % 10000 == 0:
                print(i, arr[:4])
                # break
    return stat


def conv_genejson2table(jsonfile, tabfile):
    ds = ""
    with open(jsonfile) as jfp:
        ds = json.load(jfp)
    # stat = load_filed_stat_from_datafile()

    f = open(tabfile, 'w')
    cont = []
    no = 0
    for s1 in ds['gene_source']:
        datafile = os.path.join(ds['datafile_path'], s1['datafile'])
        # stat = load_field_stat(datafile, s1['fields'], s1['format'])
        for f1 in s1['fields']:
            if not ('is_available' in f1.keys() and f1['is_available'] == False):
                no += 1
                fieldname = f1['name']
                if 'name2' in f1.keys():
                    fieldname = f1['name2']
                # fieldname = s1['name'].lower() + "_" + fieldname.lower()
                cont = [str(no), fieldname]

                cont.append(s1['name'])
                cont.append(s1['version'])
                cont.append(dv(f1, 'type'))
                is_list = 'N'
                if 'is_list' in f1.keys() and f1['is_list'] == True:
                    is_list = 'Y'
                cont.append(is_list)
                # cont.append(dv(f1, 'delimiter'))
                # k1 = fieldname
                # try:
                #     maxlen = stat['maxlen'][k1]
                # except KeyError:
                #     maxlen = 999
                # cont.append(str(maxlen))

                # values = []
                # try:
                #     values = list(stat['cnt'][k1].keys())
                # except KeyError:
                #     values = []
                cont.append(dv(f1, 'desc'))
                exvalue = ''
                # for v1 in values[:20]:
                #     if v1 != '':
                #         exvalue += ' ' + v1
                # cont.append(exvalue.strip().replace(' ',';'))
                # cont.append(';'.join(values[:20]))
                f.write('\t'.join(cont) + '\n')
    f.close()

if __name__ == "__main__":
    import proc_util
    import file_util
    import struct_util

    no_datafile = 30

    jsonfile = "/home/mk446/bio/mutanno/SRC/tests/datastructure_v0.3.0_mvp.json"
    tabfile = jsonfile + ".tab"
    datafile = "/home/mk446/bio/mutanno/DATASOURCE/MUTANOANNOT/mvp_gene_datasource_v0.3.bed"
    conv_genejson2table(jsonfile, tabfile)
