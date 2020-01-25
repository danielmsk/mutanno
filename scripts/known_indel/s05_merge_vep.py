#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s05_merge_vep.py
# made by Daniel Minseok Kwon
# 2020-01-23 11:39:39
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


def s05_merge_vep():
    flagheader = True
    tchrom = ''
    header = ''
    f = ''
    for k in range(1, 7000):
        no = str_util.zero_format(k, 6)
        inputvcf = path + title + '_' + no + '.vcf'
        print(inputvcf)
        if file_util.is_exist(inputvcf):
            vep = inputvcf + '.vep.txt'
            print("Processing..", vep)

            if not check_vcf_vep(inputvcf, vep):
                print('Error:', vep)
                exit()

            for line in open(vep):
                if line[0] == '#':
                    if flagheader:
                        if line[:len('#CHROM')] == '#CHROM':
                            arr = line.split('\t')
                            cont = arr[0]
                            cont += '\t' + arr[1]
                            cont += '\t' + arr[2]
                            cont += '\t' + arr[3]
                            cont += '\t' + arr[4]
                            cont += '\t' + arr[7].strip()
                            cont += '\n'
                            line = cont
                        header += line
                else:
                    arr = line.split('\t')
                    chrom = arr[0]
                    if chrom != tchrom:
                        try:
                            f.close()
                        except AttributeError:
                            pass
                        out2 = out.replace('#CHROM#', chrom)
                        print(out2)
                        f = open(out2, 'w')
                        tchrom = chrom
                        f.write(header)
                        flagheader = False
                    cont = arr[0].replace('chr', '')
                    cont += '\t' + arr[1]
                    cont += '\t' + arr[2]
                    cont += '\t' + arr[3]
                    cont += '\t' + arr[4]
                    cont += '\t' + arr[7].strip()
                    cont += '\n'
                    f.write(cont)
    f.close()


if __name__ == "__main__":
    import file_util
    import str_util
    title = 'sindel'
    path = '/home/mk446/mutanno/DATASOURCE/KNOWN_INDEL/hg38/tmp/'
    out = "/home/mk446/mutanno/DATASOURCE/KNOWN_INDEL/hg38/known_indel.#CHROM#.vep.tsi"
    s05_merge_vep()
