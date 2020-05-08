#!/usr/bin/env python
# -*- coding: utf-8 -*-
# stat_ensemble_gtf.py
# made by Daniel Minseok Kwon
# 2020-04-29 16:09:43
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
import test_consequence
def parse_fields(arr):
    d = {}
    d['chrom'] = arr[0].strip()
    d['spos'] = int(arr[3].strip())
    d['epos'] = int(arr[4].strip())
    d['type'] = arr[2].strip()
    d['strand'] = arr[6].strip()
    for f1 in arr[-1].split(';'):
        if f1.strip() != "":
            arr = f1.strip().split(' "')
            d[arr[0]] = arr[1].replace('"','')
    return d

def stat_ensemble_gtf(ensgene_gtf):
    
    stat = {}
    target_fields = ['type', 'gene_version', 'gene_source', 'gene_biotype', 'transcript_support_level', 'transcript_biotype']

    for tf in target_fields:
        stat[tf] = {}

    for line in file_util.gzopen(ensgene_gtf):
        line = file_util.decodeb(line)
        if line[0] != '#':
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            d = parse_fields(arr)
            
            for tf in target_fields:
                if tf in d.keys():
                    try:
                        stat[tf][d[tf]] += 1
                    except KeyError:
                        stat[tf][d[tf]] = 1
            # break
    cont = ""
    for tf in stat.keys():
        cont += tf + '\n'
        for c1 in stat[tf].keys():
            cont += '\t' + c1 +'\t'+ str(stat[tf][c1]) + '\n'
    
    file_util.fileSave(ensgene_gtf + '.stat', cont ,'w')
    print('Saved', ensgene_gtf + '.stat')

    print(cont)
            
    

if __name__ == "__main__":
    import proc_util
    import file_util
    ensgene_gtf = "/home/mk446/mutanno/DATASOURCE/ENSEMBL/hg38/Homo_sapiens.GRCh38.99.sorted.gtf.gz"
    stat_ensemble_gtf(ensgene_gtf)
