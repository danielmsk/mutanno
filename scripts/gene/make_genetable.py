#!/usr/bin/env python
# -*- coding: utf-8 -*-
# make_genetable.py
# made by Daniel Minseok Kwon
# 2020-02-24 13:07:45
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

def get_dict_value(dict1, key1, default):
    rst = default
    if key1 in dict1.keys():
        rst = dict1[key1]
    return rst


def load_example_values():
    i = 0
    exmap = {}
    for line in file_util.gzopen(datasource):
        line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        if line[0] == "#":
            header = arr
            header[0] = header[0][1:]
        else:
            flag_add = False
            for k in range(len(arr)):
                h1 = header[k]

                try:
                    cnt = len(exmap[h1])
                except KeyError:
                    cnt = 0

                if cnt < 10:
                    flag_add = True

                if arr[k] != "" and cnt < 10:
                    try:
                        tmp = exmap[h1]
                    except KeyError:
                        exmap[h1] = []
                        pass
                    
                    if arr[k] not in exmap[h1]:
                        exmap[h1].append(arr[k])

            if not flag_add:
                # break
                pass
    # print(line)
    return exmap

def make_genetable(datasource, json, out):

    exmap = load_example_values()

    reservedlist = ['chrom','spos','epos','ensgid']

    f = open(out,'w')
    cont = ["no","Field Name","Datasource","Version","Title","do import","Link","is_list","Description","value_example (separated by ';')"]
    f.write('\t'.join(cont)+'\n')
    ds = file_util.jsonOpen(json)
    no = 0
    for s1 in ds['gene_source']:
        if get_dict_value(s1, 'is_available', True):
            for f1 in s1['fields']:
                fieldname = f1['name']
                name2 = get_dict_value(f1, 'name2', '')
                if name2 != '':
                    fieldname = name2

                do_import = "Y"
                if not get_dict_value(f1, 'do_import', True):
                    do_import = "N"

                do_import = "Y"
                if not get_dict_value(f1, 'do_import', True):
                    do_import = "N"

                is_list = "N"
                if get_dict_value(f1, 'is_list', False):
                    is_list = "Y"

                if get_dict_value(f1, 'is_available', True):
                    if (fieldname != 'ensgid' and not get_dict_value(s1,'is_ensemblegene',False)) or get_dict_value(s1,'is_ensemblegene',False):
                        no += 1
                        cont = []
                        cont.append(str(no))
                        cont.append(fieldname)
                        cont.append(s1['name'])
                        cont.append(s1['version'])
                        cont.append(get_dict_value(f1, 'title', ''))
                        cont.append(do_import)
                        cont.append(get_dict_value(f1, 'link', ''))
                        cont.append(is_list)
                        cont.append(get_dict_value(f1, 'desc', ''))
                        if fieldname in reservedlist:
                            k1 = fieldname
                        else:
                            k1 = s1['name'] + '_' + fieldname
                            # k1 = s1['name'] + '_' + f1['name']
                        cont.append('; '.join(exmap[k1]))

                        # print('\t'.join(cont))
                        f.write('\t'.join(cont)+'\n')
    
    f.close()
    print('Saved', out)
    
    

if __name__ == "__main__":
    import proc_util
    import file_util
    datasource = "/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/mvp_gene_datasource_v0.3.2.bed.gz"
    json = "/home/mk446/mutanno/SRC/tests/datastructure_v0.3.2_mvp.json"
    

    datasource = sys.argv[1]
    json = sys.argv[2]
    out = json + '.genetab.txt'
    make_genetable(datasource, json, out)
