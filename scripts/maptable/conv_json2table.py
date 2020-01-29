#!/usr/bin/env python
# -*- coding: utf-8 -*-
# conv_json2table.txt
# made by Daniel Minseok Kwon
# 2020-01-29 09:15:52
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

def zero_format(i, size):
    # 00001000 format
    s1 = "{:0>" + str(size) + "d}"
    return s1.format(i)

def load_field_stat(datafile, fields, fileformat):
    stat = {}
    stat['maxlen'] = {}
    stat['cnt'] = {}
    i = 0
    h = []
    for line in file_util.gzopen(datafile.replace('#CHROM#','1')):
        line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        if line[0] == '#':
            arr[0] = arr[0][1:]
            if fileformat == 'tsi':
                h = arr[-1].split('|')
            else:
                h = arr
            # for k in range(len(arr)):
            #     h[arr[k]] = k
        else:
            if fileformat == 'tsi':
                values = arr[-1].split(',')[0].split('|')
            else:
                values = arr
            for k in range(len(values)):
                try:
                    if len(values[k]) > stat['maxlen'][h[k]]:
                        stat['maxlen'][h[k]] = len(values[k])
                except KeyError:
                    stat['maxlen'][h[k]] = len(values[k])

                try:
                    tmp = stat['cnt'][h[k]]
                except KeyError:
                    stat['cnt'][h[k]] = {}

                try:
                    stat['cnt'][h[k]][values[k]] += 1
                except KeyError:
                    stat['cnt'][h[k]][values[k]] = 1
        i += 1
        if i % 10000 == 0:
            break
    return stat

def load_filed_stat_from_datafile():
    stat = {}
    stat['maxlen'] = {}
    stat['cnt'] = {}
    h = {}
    for k in range(1,no_datafile + 1):
        k2 = k * 12
        # k2 = 1
        dfile = datafile.replace('#K#', zero_format(k2, 5))
        print(k, dfile)
        i = 0
        for line in file_util.gzopen(dfile):
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            if line[0] == "#":
                for sourcefield in arr[-1].split(';'):
                    arr2 = sourcefield.split('=')
                    sourcename = arr2[0]
                    fields = arr2[1].split('|')
                    h[sourcename] = fields
                    # print(sourcename, fields)                
                    # fkey = sourcename.lower() + "_" + f1.lower()
            else:
                i += 1
                for sourcefield in arr[-1].split(';'):
                    arr2 = sourcefield.split('=')
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

def save_stat(stat):
    for k1 in stat['cnt'].keys():
        (ks, vs) = struct_util.sortdict(stat['cnt'][k1], "desc")
        
        cntfile = './stat_tmp/' + k1 + '.stat'
        cont = 'MaxLen:' + str(stat['maxlen'][k1]) + '\n\n'
        for f1 in ks:
            cont += f1 + '\t' + str(stat['cnt'][k1][f1]) + '\n'
        file_util.fileSave(cntfile, cont, 'w')


def dv(d1, k1):
    v1 = ''
    if k1 in d1.keys():
        v1 = d1[k1]
    return v1


def conv_json2table(jsonfile, tabfile):
    ds = ""
    with open(jsonfile) as jfp:
        ds = json.load(jfp)

    stat = load_filed_stat_from_datafile()
    save_stat(stat)

    f = open(tabfile, 'w')
    cont = []
    no = 0
    for s1 in ds['source']:
        datafile = os.path.join(ds['datafile_path'], s1['datafile'])
        # stat = load_field_stat(datafile, s1['fields'], s1['format'])
        for f1 in s1['fields']:
            if not ('is_available' in f1.keys() and f1['is_available'] == False):
                no += 1
                fieldname = f1['name']
                if 'name2' in f1.keys():
                    fieldname = f1['name2']
                fieldname = s1['name'].lower() + "_" + fieldname.lower()
                cont = [str(no), fieldname]
                cont.append(s1['name'])
                cont.append(s1['version'])
                cont.append(dv(f1, 'type'))
                is_list = 'N'
                if 'is_list' in f1.keys() and f1['is_list'] == True:
                    is_list = 'Y'
                cont.append(is_list)
                # cont.append(dv(f1, 'delimiter'))
                k1 = fieldname
                try:
                    maxlen = stat['maxlen'][k1]
                except KeyError:
                    maxlen = 999
                cont.append(str(maxlen))

                values = []
                try:
                    values = list(stat['cnt'][k1].keys())
                except KeyError:
                    values = []
                cont.append(dv(f1, 'desc'))
                exvalue = ''
                for v1 in values[:20]:
                    if v1 != '':
                        exvalue += ' ' + v1
                cont.append(exvalue.strip().replace(' ',';'))
                # cont.append(';'.join(values[:20]))
                f.write('\t'.join(cont) + '\n')
    f.close()

if __name__ == "__main__":
    import proc_util
    import file_util
    import struct_util

    no_datafile = 30

    jsonfile = "/home/mk446/mutanno/SRC/tests/datastructure_v0.3.0_mvp.json"
    tabfile = jsonfile + ".tab"
    datafile = "/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/mvp_datasource_v0.3_test_tmp/#K#.tsi"
    conv_json2table(jsonfile, tabfile)
    
    # print(stat)
