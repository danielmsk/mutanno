#!/usr/bin/env python
# -*- coding: utf-8 -*-
# /home/mk446/bio/mutanno/SRC/scripts/gene/preproc_GTEx.py
# made by Daniel Minseok Kwon
# 2020-04-20 14:52:03
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
import time


def preproc_GTEx(gtex_file, out):
    f = open(out, 'w')
    for line in file_util.gzopen(gtex_file):
        line = file_util.decodeb(line)
        arr = line.split('\t')
        cont = []
        if arr[0] == "#Name":
            cont = ['#ensgid', 'gtex_expression']
        elif not 'PAR_Y' in arr[0]:
            ensgid = arr[0].split('.')[0]
            cont = [ensgid, arr[0].strip()]
        if len(cont) == 2:
            f.write('\t'.join(cont) + '\n')
        # break
    f.close()

    time.sleep(3)
    proc_util.run_cmd('gz ' + out)


if __name__ == "__main__":
    import proc_util
    import file_util
    url = "https://gtexportal.org/home/gene/<ID>"
    gtex_file = "../../../DATASOURCE/EXPRESSION/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.mod.gz"
    out = "../../../DATASOURCE/EXPRESSION/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_link.tsv"
    preproc_GTEx(gtex_file, out)
