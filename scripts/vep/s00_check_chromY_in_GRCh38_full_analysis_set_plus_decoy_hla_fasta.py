#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s00_check_chromY_in_GRCh38_full_analysis_set_plus_decoy_hla_fasta.py
# made by Daniel Minseok Kwon
# 2020-02-27 15:39:08
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


def s00_check_chromY_in_GRCh38_full_analysis_set_plus_decoy_hla_fasta():
    flag = False
    for line in open(fasta):
        if line[0] == '>':
            if flag:
                flag = False

            if line[:len('>chrY')] == '>chrY':
                flag = True

            if line[:len('>chrM')] == '>chrM':
                flag = False

            if line[:len('>chr1_KI270706v1_random')] == '>chr1_KI270706v1_random':
                flag = False



        if flag:
            print(line.strip())

            # break

def print_base(chrfa, p1):
    cont = ''
    for line in open(chrfa):
        if line[0] != '>':
            cont += line.strip()
    print(cont[p1])

if __name__ == "__main__":
    import proc_util
    import file_util
    fasta = "/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh38d1/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    # fasta = "/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/hg38.ucsc/ucsc.hg38.sorted.fa"
    s00_check_chromY_in_GRCh38_full_analysis_set_plus_decoy_hla_fasta()
    # print_base('chrY_ucsc.fa', 57000001)