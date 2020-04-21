#!/usr/bin/env python
# -*- coding: utf-8 -*-
# preproc_BrainAtlas_link.py
# made by Daniel Minseok Kwon
# 2020-04-20 15:27:11
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
sys.path.append("..")


def preproc_BrainAtlas_link(out, probe_file):
    genesymbol2ensgid, ensgid2genesymbol = preproc_util.get_map_genesymbol_ensgid()

    f = open(out, 'w')
    for line in file_util.gzopen(probe_file):
        line = file_util.decodeb(line)
        arr = line.split('\t')
        gene_symbol = arr[3].strip()
        if line[0] == "#":
            ensgid = "#ensgid"
            gene_symbol = "brainatlas_microarray"
        else:
            try:
                ensgid = genesymbol2ensgid[gene_symbol]
            except KeyError:
                ensgid = ""

        if ensgid != "":
            cont = [ensgid, gene_symbol]
            f.write('\t'.join(cont) + '\n')

    f.close()
    time.sleep(3)
    proc_util.run_cmd('gz ' + out, True)
    

if __name__ == "__main__":
    import proc_util
    import file_util
    import preproc_util
    path = preproc_util.DATASOURCEPATH + "/EXPRESSION/BRAINATLAS/"
    probe_file = path + "H0351.1009.microarray/Probes.tsv.gz"
    out = path + "BRAINATLAS_microarry_link.tsv"
    preproc_BrainAtlas_link(out, probe_file)
