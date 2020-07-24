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

def gv(d1, k1):
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


fieldlist = []
fieldlist.append('ensgid')
fieldlist.append('gene_version')
fieldlist.append('gene_name')
fieldlist.append('gene_source')
fieldlist.append('gene_biotype')
fieldlist.append('transcript_id')
fieldlist.append('transcript_version')
fieldlist.append('transcript_name')
fieldlist.append('transcript_source')
fieldlist.append('transcript_biotype')
fieldlist.append('transcript_support_level')
fieldlist.append('tag')
fieldlist.append('exon_id')
fieldlist.append('exon_number')
fieldlist.append('exon_version')
fieldlist.append('protein_id')
fieldlist.append('protein_version')
fieldlist.append('ccds_id')


def convert_ensemblgene_gtf2bed(gtf, bed, gtftypes, filetype="bed", biotypes=[]):
    global fieldlist

    print(gtf)

    f = open(bed, 'w')

    if filetype == "bed":
        cont = ["#CHROM", "SPOS", "EPOS", "type" ,"strand"]
        cont.extend(fieldlist)
    elif filetype == "tsi":
        cont = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
        info = "ENSEMBLANNOT=type|strand|" + '|'.join(fieldlist)
        cont.append(info)
    f.write('\t'.join(cont) + '\n')
    i = 0
    cnt = {}

    # print(len(fieldlist))
    for line in file_util.gzopen(gtf):
        line = line.decode('UTF-8')
        if line[0] != '#':
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()

            gtftype = arr[2]

            if (gtftypes[0] == "alltype") or (gtftype in gtftypes):
                i += 1
                if arr[0] == "MT":
                    arr[0] = "M"
                
                cnt = count_no_fields(cnt, 'chrom', arr[0])
                m = {}
                for f1 in arr[-1].strip().split(';'):
                    arr2 = f1.strip().split(' ')
                    fieldname = arr2[0].strip()
                    value = ' '.join(arr2[1:]).replace('"','')
                    if fieldname == "gene_id":
                        fieldname = "ensgid"
                    m[fieldname] = value

                    if fieldname not in ['gene_id','gene_name']:
                        cnt = count_no_fields(cnt, fieldname, value)

                if filetype == "bed":
                    if (len(biotypes) == 0) or (m['gene_biotype'] in biotypes):
                        cont = [arr[0], str(int(arr[3])-1), arr[4], arr[2] , arr[6]]
                        for f1 in fieldlist:
                            cont.append(gv(m,f1))
                        f.write('\t'.join(cont) + '\n')
                elif filetype == "tsi":
                    cont = [arr[0], arr[3], "", "", "", "", ""]
                    
                    info2 = ""
                    for f1 in fieldlist:
                        info2 += "|" + gv(m,f1)

                    # if info2.replace("|","").strip() != "":
                    if (len(biotypes) == 0) or (m['gene_biotype'] in biotypes):
                        info = "END=" + arr[4]
                        info += ";ENSEMBLANNOT=" + gtftype + "|"
                        if arr[6] == "+":
                            info += "1"
                        else:
                            info += "0"
                        info += info2
                        cont.append(info)
                        # print(len(info.split('|')))
                        f.write('\t'.join(cont) + '\n')
                        # break
                if i % 10000 == 0:
                    print(i, arr[:4])
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
    # convert_ensemblgene_gtf2bed(gtf, path + "Homo_sapiens.GRCh38.99.bed", ["gene"], "bed")
    # convert_ensemblgene_gtf2bed(gtf, path + "Homo_sapiens.GRCh38.99.alltype.bed", ["alltype"], "bed")
    convert_ensemblgene_gtf2bed(gtf, path + "Homo_sapiens.GRCh38.99.codeing_mirna_polymorphic.bed", ["alltype"], "bed", ['protein_coding','miRNA','polymorphic_pseudogene'])
    # convert_ensemblgene_gtf2bed(gtf, path + "Homo_sapiens.GRCh38.99.codeing_mirna_polymorphic.tsi", ["alltype"], "tsi", ['protein_coding','miRNA','polymorphic_pseudogene'])
