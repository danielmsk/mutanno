#!/usr/bin/env python
# -*- coding: utf-8 -*-
# preproc_BRAINATLAS.py
# made by Daniel Minseok Kwon
# 2020-03-24 14:05:42
#########################
import preproc_util
import sys
import os
import csv
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path = "/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path = "/ms1/bin/python_lib"
else:
    sys_path = "/home/mk446/bin/python_lib"
sys.path.append(sys_path)

sys.path.append("..")


def load_ontology(csvfile):
    m = {}
    h = {}
    with open(csvfile, newline='') as cf:
        for arr in csv.reader(cf, delimiter=',', quotechar='"'):
            if arr[0] == "id":
                for k in range(len(arr)):
                    h[arr[k]] = k
            else:
                m[arr[h['id']]] = arr[h['name']]
    return m


def load_sampleannot(csvfile, ontology):
    m = []
    h = {}
    for line in open(csvfile):
        arr = line.split(',')
        arr[-1] = arr[-1].strip()
        if arr[0] == "RNAseq_sample_name":
            for k in range(len(arr)):
                h[arr[k]] = k
        else:
            ontology_name = ontology[arr[h['ontology_structure_id']]]
            m.append(ontology_name)
            # print (arr[h['main_structure']], arr[h['sub_structure']], ontology_name)
        # print(arr)
    return m


def preproc_BRAINATLAS(dpath, tpm_csv, annot_csv, ontology_csv, out):
    f = open(dpath + out, 'w')
    ontology = load_ontology(dpath + ontology_csv)
    sampleannot = load_sampleannot(dpath + annot_csv, ontology)

    cont = []
    cont.append('#ENSG')
    cont.append('GeneSymbol')
    cont.extend(sampleannot)
    f.write('\t'.join(cont) + '\n')

    map_gene_ensgid = preproc_util.get_map_genesymbol_ensgid_from_hgnc()
    for line in open(dpath + tpm_csv):
        arr = line.split(',')
        arr[-1] = arr[-1].strip()
        genesymbol = arr[0].replace('"', '')
        try:
            ensgid = map_gene_ensgid[genesymbol]
            # print(genesymbol, ensgid, arr[1:])
            cont = [ensgid]
            cont.append(genesymbol)
            cont.extend(arr[1:])
            f.write('\t'.join(cont) + '\n')
        except KeyError:
            pass
            # print(genesymbol)
        # break
    # print(map_gene_ensid)
    f.close()
    print('Saved', dpath + out)


if __name__ == "__main__":
    path = preproc_util.DATASOURCEPATH + "/EXPRESSION/BRAIN_ATLAS/"
    out = "RNAseqTPM.tsv"
    # preproc_BRAINATLAS(path + "rnaseq_donor10021/", "RNAseqTPM.csv", "SampleAnnot.csv", "Ontology.csv", out)
    preproc_BRAINATLAS(path + "rnaseq_donor9861/", "RNAseqTPM.csv", "SampleAnnot.csv", "Ontology.csv", out)
