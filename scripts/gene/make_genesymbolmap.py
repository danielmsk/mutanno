#!/usr/bin/env python
# -*- coding: utf-8 -*-
# make_genesymbolmap.py
# made by Daniel Minseok Kwon
# 2020-03-31 14:24:31
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
sys.path.append('..')


def make_genesymbolmap(out):
    f = open(out, 'w')
    genesymbolmap = preproc_util.get_map_genesymbol_ensgid_from_ensembl(flag_upper=False)
    genesymbolmap2 = preproc_util.get_map_genesymbol_ensgid_from_hgnc(flag_upper=False)
    
    ensgmap = {}
    for gsymbol in genesymbolmap.keys():
        ensgid = genesymbolmap[gsymbol]
        if ensgid not in ensgmap.keys():
            ensgmap[ensgid] = [gsymbol]
        else:
            print('ENSML',ensgid, gsymbol, ensgmap[ensgid])
            ensgmap[ensgid].append(gsymbol)

    for gsymbol in genesymbolmap2.keys():
        ensgid = genesymbolmap2[gsymbol]
        if ensgid not in ensgmap.keys():
            ensgmap[ensgid] = [gsymbol]
        else:
            # print('HGNC', ensgid, gsymbol, ensgmap[ensgid])
            ensgmap[ensgid].append(gsymbol)

    h = {}
    for line in file_util.gzopen(preproc_util.NCBIgene):
        line = file_util.decodeb(line)
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        if line[0] == "#":
            arr[0] = arr[0][1:]
            for k in range(len(arr)):
                h[arr[k]] = k
        else:
            ensgid = ""
            for f1 in arr[h['dbXrefs']].split('|'):
                if 'Ensembl:' in f1:
                    ensgid = f1.split(':')[-1].strip()

            if ensgid != "":
                symlist = [arr[h['Symbol']].strip()]
                symlist.extend(arr[h['Synonyms']].strip().split('|'))
                for gsymbol in symlist:
                    if gsymbol != "" and gsymbol != "-":
                        if ensgid not in ensgmap.keys():
                            ensgmap[ensgid] = [gsymbol]
                        else:
                            # print('HGNC', ensgid, gsymbol, ensgmap[ensgid])
                            ensgmap[ensgid].append(gsymbol)


    for ensgid in ensgmap.keys():
        cont = []
        cont.append(ensgid)
        cont.append(ensgmap[ensgid][0])
        cont.append('|'.join(ensgmap[ensgid]))
        f.write('\t'.join(cont) + '\n')
    f.close()
    

if __name__ == "__main__":
    import preproc_util
    import proc_util
    import file_util
    path = preproc_util.DATASOURCEPATH + "/GENE/"
    out = preproc_util.DATASOURCEPATH + "/GENE/gene_symbolmap.tsv" 
    make_genesymbolmap(out)
