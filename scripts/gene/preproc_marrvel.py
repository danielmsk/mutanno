#!/usr/bin/env python
# -*- coding: utf-8 -*-
# preproc_marrvel.py
# made by Daniel Minseok Kwon
# 2020-03-31 04:55:56
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


def preproc_marrvel(out):
    f = open(out, 'w')
    ensg_map = preproc_util.get_map_genesymbol_ensgid_from_ensembl()
    cont = ['#ensgid','Marrvel']
    f.write('\t'.join(cont) + '\n')
    for k1 in ensg_map.keys():
        cont = [ensg_map[k1], k1]
        f.write('\t'.join(cont) + '\n')
    f.close()
    print('Saved',out)

    

if __name__ == "__main__":
    import preproc_util
    path = preproc_util.DATASOURCEPATH + "/GENE/MARRVEL/"
    out = path + "marrvel_genesymbol.tsv"
    preproc_marrvel(out)
