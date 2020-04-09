#!/usr/bin/env python
# -*- coding: utf-8 -*-
# convert_ensemblgene_gtf2bed.py
# made by Daniel Minseok Kwon
# 2020-02-03 16:40:07
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

def get_value_from_dict(d1, k1):
    v1 = ''
    try:
        v1 = d1[k1]
    except KeyError:
        pass
    return v1

def count_no_fields(cnt, fieldname, value):
    try:
        tmp = cnt[fieldname]
    except KeyError:
        cnt[fieldname] = {}

    try:
        cnt[fieldname][value] += 1
    except KeyError:
        cnt[fieldname][value] = 1
    return cnt

def convert_ensemblgene_gtf2bed(gtf, bed):
    f = open(bed, 'w')
    cont = ["#chrom", "spos", "epos", "strand", "ensgid", "gene_version", "gene_symbol"]
    cont.extend(["gene_source", "gene_biotype"])
    f.write('\t'.join(cont) + '\n')
    i = 0
    cnt = {}
    for line in file_util.gzopen(gtf):
        line = line.decode('UTF-8')
        if line[0] != '#':
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            if arr[2] == "gene":
                i += 1
                if arr[0] == "MT":
                    arr[0] = "M"
                cont = [arr[0], arr[3], arr[4], arr[6]]
                cnt = count_no_fields(cnt, 'chrom', arr[0])
                m = {}
                for f1 in arr[-1].strip().split(';'):
                    arr2 = f1.strip().split(' ')
                    fieldname = arr2[0].strip()
                    value = ' '.join(arr2[1:]).replace('"','')
                    m[fieldname] = value

                    if fieldname not in ['gene_id','gene_name']:
                        cnt = count_no_fields(cnt, fieldname, value)
                cont.append(get_value_from_dict(m,'gene_id'))
                cont.append(get_value_from_dict(m,'gene_version'))
                cont.append(get_value_from_dict(m,'gene_name'))
                cont.append(get_value_from_dict(m,'gene_source'))
                cont.append(get_value_from_dict(m,'gene_biotype'))
                # cont = arr
                f.write('\t'.join(cont) + '\n')
                if i > 100:
                    # break
                    pass
    f.close()

    cont = ''
    for k1 in cnt.keys():
        # print(k1)
        cont += k1 + '\n'
        for k2 in cnt[k1].keys():
            cont += '\t' + k2 + '\t' + str(cnt[k1][k2]) + '\n'
            # print('\t',k2, cnt[k1][k2])
    file_util.fileSave(bed + '.stat', cont, 'w')
    proc_util.run_cmd('tabixgzbed ' + bed, True)

if __name__ == "__main__":
    import proc_util
    import file_util
    path = "/home/mk446/bio/mutanno/DATASOURCE/ENSEMBL/hg38/"
    gtf = path + "Homo_sapiens.GRCh38.99.gtf.gz"
    bed = path + "Homo_sapiens.GRCh38.99.bed"
    convert_ensemblgene_gtf2bed(gtf, bed)
