#!/usr/bin/env python
# -*- coding: utf-8 -*-
# preproc_dbNSFP.py
# made by Daniel Minseok Kwon
# 2020-03-05 08:54:00
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

def change_dot_to_blank(arr):
    for k in range(len(arr)):
        arr[k] = arr[k].strip()
        if arr[k] == ".":
            arr[k] = ""
    return arr

def cleanup_last_empty_field(arr):
    for k in range(len(arr)):
        if arr[k] != '' and arr[k][-1] == ';':
            arr[k] = arr[k][:-1]
    return arr


def split_gwas_trait_id_and_desc(trait_desc):
    trait_list = []
    pmid_list = []
    if trait_desc != "":
        for t1 in trait_desc.split('];'):
            arr2 = t1.split('[')
            # print(trait_desc, arr2)
            trait_list.append(arr2[0].strip())
            pmid_list.append(arr2[1].replace(']','').replace(';','~'))
    return ';'.join(trait_list), ';'.join(pmid_list)

def preproc_dbNSFP(infile, outfile):
    f = open(outfile, 'w')
    i = 0
    h = {}
    for line in file_util.gzopen(infile):
        line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        if i == 0:
            for k in range(len(arr)):
                h[arr[k]] = k
            arr.append('Trait_association(GWAS)_PMID')
            f.write('#'+'\t'.join(arr) + '\n')
        else:
            arr = change_dot_to_blank(arr)
            arr = cleanup_last_empty_field(arr)
            idx = h['Trait_association(GWAS)']
            arr[idx], gwas_pmid = split_gwas_trait_id_and_desc(arr[idx])
            arr.append(gwas_pmid)
            f.write('\t'.join(arr) + '\n')
        i += 1

    f.close()

    proc_util.run_cmd('tabixgzbed ' + outfile)
    

if __name__ == "__main__":
    import proc_util
    import file_util
    path = "/home/mk446/mutanno/DATASOURCE/VARIANTDB/dbNSFP/hg38/dbNSFP4.0c/"
    infile = path + "dbNSFP4.0_gene.complete.gz"
    outfile = path + "dbNSFP4.0_gene.complete.mod"
    preproc_dbNSFP(infile, outfile)
