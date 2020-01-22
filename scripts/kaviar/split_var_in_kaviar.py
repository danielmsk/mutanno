#!/usr/bin/env python
# -*- coding: utf-8 -*-
# split_var_in_kaviar.py
# made by Min-Seok Kwon
# 2019-11-04 17:18:42
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


def split_var_in_kaviar(orifile):
    out = orifile.replace('.vcf.gz', '') + '.splited.vcf'
    f = open(out, 'w')
    for line in file_util.gzopen(orifile):
        line = line.decode('UTF-8')
        if line[0] != "#":
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            ori = arr[1] + '_' + arr[3] + '_' + arr[4]
            if ',' in arr[4]:
                arr4 = arr[4].strip().split(',')
                endf = ''
                for f1 in arr[7].split(';'):
                    arrf = f1.split('=')
                    if arrf[0] == "AC":
                        arr_ac = arrf[1].split(',')
                    if arrf[0] == "AF":
                        arr_af = arrf[1].split(',')
                    if arrf[0] == "AN":
                        an = arrf[1]
                    if arrf[0] == "END":
                        endf = arrf[1]
                for k in range(len(arr4)):
                    ref = arr[3]
                    alt = arr4[k]
                    j = 0
                    i = 0
                    p = 0
                    if min(len(ref), len(alt)) > 1:
                        for j in range(min(len(ref), len(alt))):
                            if len(ref) == 1 or len(alt) == 1:
                                break
                            if alt[0] == ref[0]:
                                alt = alt[1:]
                                ref = ref[1:]
                                p += 1
                            else:
                                break

                        for i in range(1, min(len(ref), len(alt))):
                            if len(ref) == 1 or len(alt) == 1:
                                break
                            if alt[-1] == ref[-1]:
                                alt = alt[:-1]
                                ref = ref[:-1]
                            else:
                                break
                        i += 1
                    cont = []
                    cont.append(arr[0])
                    cont.append(str(int(arr[1]) + p))
                    cont.append(arr[2])
                    cont.append(ref)
                    cont.append(alt)
                    # cont.append(ori)
                    cont.extend(arr[5:7])
                    cinfo = "AF=" + arr_af[k] + ';AC=' + arr_ac[k] + ';AN=' + an
                    if endf != '':
                        cinfo += ';END=' + endf
                    cont.append(cinfo)

                    f.write('\t'.join(cont) + '\n')
                # break
            else:
                f.write('\t'.join(arr) + '\n')
        else:
            f.write(line)

    f.close()
    proc_util.run_cmd('tabixgz ' + out)


if __name__ == "__main__":
    import file_util
    import proc_util
    path = "/home/mk446/mutanno/DATASOURCE/POPULATION_AF/KAVIAR/hg38/"
    split_var_in_kaviar(path + "Kaviar-160204-Public-hg38-trim.vcf.gz")
    # split_var_in_kaviar("test.vcf.gz")
