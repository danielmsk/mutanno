#!/usr/bin/env python
# -*- coding: utf-8 -*-
# preproc_MedGen.py
# made by Daniel Minseok Kwon
# 2020-03-06 17:03:40
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


def load_mim2gene_medgen():
    mimmap = {}
    genemap = {}
    h = {}
    for line in open(mim2gene_medgen):
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        if line[0] == '#':
            arr[0] = arr[0].replace('#','')
            for k in range(len(arr)):
                h[arr[k]] = k
        else:
            mimmap[arr[h['MIM number']]] = arr[h['MedGenCUI']].strip()
            genemap[arr[h['GeneID']]] = arr[h['MedGenCUI']].strip()
    return mimmap, genemap

def load_omim_medgene():
    mimmap = {}
    genemap = {}
    h = {}
    for line in file_util.gzopen(omim_medgene_file):
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        if line[0] == '#':
            arr[0] = arr[0].replace('#','')
            for k in range(len(arr)):
                h[arr[k]] = k
        else:
            mimmap[arr[h['MIM number']]] = arr[h['MedGenCUI']].strip()
            genemap[arr[h['GeneID']]] = arr[h['MedGenCUI']].strip()
    return mimmap, genemap


def preproc_MedGen(out):
    # mimmap, genemap = load_mim2gene_medgen()
    mimmap, genemap = load_omim_medgene()

    out = ''
    h = {} 
    for line in file_util.gzopen(dbNSFP_file):
        line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        if line[0] == '#':
            arr[0] = arr[0].replace('#','')
            for k in range(len(arr)):
                h[arr[k]] = k
        else:
            # arr[h['MIM_phenotype_id']]
            ensg_id = arr[h['Ensembl_gene']]
            mim_phenoids = arr[h['MIM_phenotype_id']]
            if mim_phenoids != "":
                for mim_phenoid in mim_phenoids.split(';'):
                    print(mim_phenoid)
                    print(mimmap[mim_phenoid])


if __name__ == "__main__":
    import proc_util
    import file_util
    path = "/home/mk446/mutanno/DATASOURCE/NCBI/MEDGEN/"
    # mim2gene_medgen = path + "mim2gene_medgen"

    omim_medgene_file = path + "MedGen_HPO_OMIM_Mapping.txt.gz"
    out = path + "MedGen_CUI_ENSG.tsv"
    #MIM_phenotype_id
    dbNSFP_file = "/home/mk446/mutanno/DATASOURCE/VARIANTDB/dbNSFP/hg38/dbNSFP4.0c/dbNSFP4.0_gene.complete.mod.gz"
    preproc_MedGen(out)
