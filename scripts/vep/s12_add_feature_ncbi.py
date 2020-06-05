#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s12_add_feature_ncbi.py
# made by Daniel Minseok Kwon
# 2020-06-01 15:38:42
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

entrezmap = {}
refseqmap = {}

def load_entrez_refseq_id():
    path = "/home/mk446/bio/mutanno/DATASOURCE/ENSEMBL/hg38/"
    entrezmap = {}
    for line in file_util.gzopen(path + 'Homo_sapiens.GRCh38.99.entrez.tsv.gz'):
        line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        entrezmap[arr[0].strip()] = arr[3].strip()
    refseqmap = {}
    for line in file_util.gzopen(path + 'Homo_sapiens.GRCh38.99.refseq.tsv.gz'):
        line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        enst_id = arr[1].strip()
        refseqmap[enst_id] = arr[3].strip()
    return entrezmap, refseqmap


def get_refseq_id(transcriptid):
    global entrezmap, refseqmap
    if len(refseqmap.keys()) == 0:
        entrezmap, refseqmap = load_entrez_refseq_id()
    try:
        tid = refseqmap[transcriptid]
    except KeyError:
        tid = ''
    return tid


def s12_add_feature_ncbi(chrom):
    vepfile = "/home/mk446/mutanno/DATASOURCE/ANNOT/VEP/hg38/v99/vep.99.hg38."+chrom+".sorted.rmcsq.tsi.gz"
    out = vepfile + '.mod'
    f = open(out, 'w')
    i = 0
    for line in file_util.gzopen(vepfile):
        line = file_util.decodeb(line)
        i += 1
        if line[:2] != "##":
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            if line[0] == "#":
                arr[-1] += "|Feature_ncbi"
                col = arr[-1].split('|')
                colnum1 = len(col)
            else:
                sections = []
                for sec in arr[-1].split(','):
                    arr2 = sec.split('|')
                    sec += "|" + get_refseq_id(arr2[col.index('Feature')])
                    colnum2 = len(sec.split('|'))
                    if colnum1 != colnum2:
                        print(line)
                        raise Exception("Not matched column number.")
                    sections.append(sec)
                arr[-1] = ','.join(sections)
            line = '\t'.join(arr) + '\n'

        f.write(line)
        if i % 4000000 == 0:
            print(i, arr[:5], colnum1, colnum2)
    f.close()


def run():
    for chrom in seq_util.MAIN_CHROM_LIST:
        cmd = "python s12_add_feature_ncbi.py " + chrom
        print(cmd)

if __name__ == "__main__":
    import seq_util
    import file_util
    s12_add_feature_ncbi(sys.argv[1])
    # run()
