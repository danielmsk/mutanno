#!/usr/bin/env python
# -*- coding: utf-8 -*-
# preproc_shet.py
# made by Daniel Minseok Kwon
# 2020-03-08 21:59:44
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
sys.path.append("..")


def preproc_shet(shet_file, out):
    hgnc_ensg_map = preproc_util.get_map_genesymbol_ensgid_from_hgnc()
    hgnc_ensg_map2 = preproc_util.get_map_genesymbol_ensgid_from_ensembl()

    # print(hgnc_ensg_map)

    f = open(out, 'w')

    for line in open(shet_file):
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        if arr[0] == "gene_symbol":
            cont = ['#ensgid']
            cont.extend(arr)
        else:
            gene_symbol = arr[0].strip()
            try:
                ensg_id = hgnc_ensg_map[gene_symbol.upper()]
            except KeyError:
                try:
                    ensg_id = hgnc_ensg_map2[gene_symbol.upper()]
                except KeyError:
                    print('KeyError: ' + gene_symbol)
                    ensg_id = ''
            cont = [ensg_id]
            for k in range(len(arr)):
                arr[k] = arr[k].strip()
            cont.extend(arr)
        f.write('\t'.join(cont) + '\n')
    f.close()
    print('Saved', out)


if __name__ == "__main__":
    import preproc_util

    path = preproc_util.DATASOURCEPATH + "/GENE/SHET/"
    preproc_shet(path + "075523-1.txt", path + "075523-1.tsv")
