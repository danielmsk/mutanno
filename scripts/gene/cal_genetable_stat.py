#!/usr/bin/env python
# -*- coding: utf-8 -*-
# make_genetable_stat.py
# made by Daniel Minseok Kwon
# 2020-02-24 13:54:43
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


def make_genetable_stat(datasource):

    reservedlist = ['spos','epos','engsid']

    cntmap = {}
    cntvalues = {}
    out = datasource + '.stat'
    f = open(out, 'w')
    i = 0
    for line in file_util.gzopen(datasource):
        line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        if line[0] == "#":
            header = arr
            header[0] = header[0][1:]
            for k in range(len(header)):
                h1 = header[k]
                cntmap[h1] = {'filled':0, 'blank':0}
        else:
            i += 1
            for k in range(len(arr)):
                h1 = header[k]
                if arr[k] == "":
                    cntmap[h1]['blank'] += 1
                else:
                    cntmap[h1]['filled'] += 1
                
                try:
                    cnt = cntmap[h1]
                except KeyError:
                    cnt = 0
            if i % 1000 == 0:
                print(i)
                # break

    cont = "Field\tFilled\tBlank"
    f.write(cont + '\n')
    for h1 in cntmap.keys():
        cont = h1 + '\t' + str(cntmap[h1]['filled']) + '\t' + str(cntmap[h1]['blank'])
        f.write(cont + '\n')
    f.close()
    print('Saved', out)
if __name__ == "__main__":
    import proc_util
    import file_util
    # datasource = "/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/mvp_gene_datasource_v0.3.1.bed.gz"
    datasource = sys.argv[1]
    make_genetable_stat(datasource)
