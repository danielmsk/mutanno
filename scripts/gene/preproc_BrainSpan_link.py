#!/usr/bin/env python
# -*- coding: utf-8 -*-
# preproc_BrainSpan_link.py
# made by Daniel Minseok Kwon
# 2020-04-20 16:21:05
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
import csv
sys.path.append("..")

def load_probe_file(genemap, probe_file, k1, tidx):
    with open(probe_file, newline='') as csvfile:
        for row in csv.reader(csvfile, delimiter=',', quotechar='"'):
            ensembl_gene_id = row[2]
            if ensembl_gene_id != "ensembl_gene_id":
                gene_symbol = row[3]
                try:
                    genemap[ensembl_gene_id][k1] = row[tidx]
                except KeyError:
                    genemap[ensembl_gene_id] = {}
                    genemap[ensembl_gene_id][k1] = row[tidx]
    return genemap


def preproc_BrainSpan_link(out, ma_probe_file, rnaseq_probe_file):

    genemap = {}
    genemap = load_probe_file(genemap, ma_probe_file, "ma", 3)
    genemap = load_probe_file(genemap, rnaseq_probe_file, "rs", 2)

    
    f = open(out, 'w')
    cont = ['#ensgid', 'brainspan_microarray', 'brainspan_rnaseq']
    f.write('\t'.join(cont) + '\n')
    for ensgid in genemap.keys():
        try:
            ma = genemap[ensgid]['ma']
        except KeyError:
            ma = ""
        try:
            rs = genemap[ensgid]['rs']
        except KeyError:
            rs = ""
        cont = [ensgid, ma, rs]
        f.write('\t'.join(cont) + '\n')

    f.close()
    time.sleep(3)
    proc_util.run_cmd('gz ' + out, True)
    

if __name__ == "__main__":
    import proc_util
    import file_util
    import preproc_util
    path = preproc_util.DATASOURCEPATH + "/EXPRESSION/BRAINSPAN/"
    ma_probe_file = path + "Developmental_Transcriptome_Dataset/Exon_microarray_summarized_to_genes_gene_array_matrix/rows_metadata.csv"
    rnaseq_probe_file = path + "Developmental_Transcriptome_Dataset/RNA-Seq_Gencode_v10_summarized_to_genes_genes_matrix/rows_metadata.csv"
    out = path + "BRAINSPAN_link.tsv"
    preproc_BrainSpan_link(out, ma_probe_file, rnaseq_probe_file)
