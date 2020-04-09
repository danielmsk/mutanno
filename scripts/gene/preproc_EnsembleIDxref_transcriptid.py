#!/usr/bin/env python
# -*- coding: utf-8 -*-
# preproc_preproc_EnsembleIDxref_transcriptid.py
# made by Daniel Minseok Kwon
# 2020-04-06 15:57:45
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


def preproc_preproc_EnsembleIDxref_transcriptid(ensxref_file, out):
    f = open(out, 'w')
    h = {}
    for line in file_util.gzopen(ensxref_file):
        line = file_util.decodeb(line)
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()

        if line[0] == "#":
            arr[0] = arr[0][1:]
            for k in range(len(arr)):
                h[arr[k]] = k
            ensgid = arr[h['EnsemblID']]
            trs = arr[h['Ensembl_TRS']]
            pro = arr[h['Ensembl_PRO']]
            cont = [ensgid, trs, pro]
            f.write('#' + '\t'.join(cont) + '\n')
        else:
            ensgid = arr[h['EnsemblID']]
            trslist = arr[h['Ensembl_TRS']].split('|')
            prolist = arr[h['Ensembl_PRO']].split('|')
            for k in range(len(trslist)):
                cont = [ensgid, trslist[k], prolist[k]]
                f.write('\t'.join(cont) + '\n')
    f.close()
    print('Saved', out)
    

if __name__ == "__main__":
    import proc_util
    import file_util
    ensxref_file = "/home/mk446/mutanno/DATASOURCE/ENSEMBL/hg38/idmap_ensembl_uniprot_xref.tsv.gz"
    out = ensxref_file.replace('.tsv.gz','') + '.transcriptid.tsv'
    preproc_preproc_EnsembleIDxref_transcriptid(ensxref_file, out)
