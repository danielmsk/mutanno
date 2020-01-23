#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s04_check_vep.py
# made by Daniel Minseok Kwon
# 2020-01-23 09:26:06
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


def check_vcf_vep(vcf, vep):
    m = {}
    for line in open(vep):
        if line[0] != '#':
            arr = line.split('\t')
            k1 = '_'.join([arr[0], arr[1], arr[3], arr[4]])
            m[k1] = 1

    flag = True
    for line in open(vcf):
        if line[0] != '#':
            arr = line.split('\t')
            k1 = '_'.join([arr[0], arr[1], arr[3], arr[4]])
            try:
                tmp = m[k1]
            except KeyError:
                print(arr)
                flag = False
    return flag


def s04_check_vep():
    i = 0
    for vcf in file_util.walk(path, '.vcf'):
        if '.vep.vcf' not in vcf:
            vep = vcf + ".vep.vcf"
            print(vep)
            if file_util.is_exist(vep):
                if not check_vcf_vep(vcf, vep):
                    print('Error:', vep)
            i += 1
    print(i)


if __name__ == "__main__":
    path = '/home/mk446/mutanno/DATASOURCE/KNOWN_INDEL/hg38/tmp3/'
    import proc_util
    import file_util
    s04_check_vep()
