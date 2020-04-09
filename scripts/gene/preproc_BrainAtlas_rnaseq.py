#!/usr/bin/env python
# -*- coding: utf-8 -*-
# preproc_BRAINATLAS.py
# made by Daniel Minseok Kwon
# 2020-03-24 14:05:42
#########################
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

def load_ontology_tsv(tsvfile):
    m = {}
    for line in open(tsvfile):
        if line[0] != '#':
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            m[arr[0]] = arr[1]
    return m

def load_sampleannot(csvfile, ontology):
    m = []
    h = {}
    acronym_list = []
    no_map = {}
    for line in file_util.gzopen(csvfile):
        line = file_util.decodeb(line)
        arr = line.split(',')
        arr[-1] = arr[-1].strip()
        print(arr)
        if arr[0] == "RNAseq_sample_name":
            for k in range(len(arr)):
                h[arr[k]] = k
        else:
            acronym = ontology[arr[h['ontology_structure_id']]]
            if acronym in acronym_list:
                try:
                    no_map[acronym] += 1
                except KeyError:
                    no_map[acronym] = 2

                # m.append(acronym + '.' + str(no_map[acronym]))
                m.append(acronym)
            else:
                m.append(acronym)
            acronym_list.append(acronym)
            # print (arr[h['main_structure']], arr[h['sub_structure']], ontology_name)
        # print(arr)
    return m


def preproc_BRAINATLAS(dpath, tpm_csv, annot_csv, ontology_tsv, out):
    f = open(dpath + out, 'w')
    ontology = load_ontology_tsv(ontology_tsv)
    sampleannot = load_sampleannot(dpath + annot_csv, ontology)

    cont = []
    cont.append('#ensgid')
    cont.append('GeneSymbol')
    cont.extend(sampleannot)
    f.write('\t'.join(cont) + '\n')

    map_gene_ensgid = preproc_util.get_map_genesymbol_ensgid_from_hgnc()
    map_gene_ensgid2 = preproc_util.get_map_genesymbol_ensgid_from_ensembl()
    for line in open(dpath + tpm_csv):
        arr = line.split(',')
        arr[-1] = arr[-1].strip()
        genesymbol = arr[0].replace('"', '')
        try:
            ensgid = map_gene_ensgid[genesymbol.upper()]
            # print(genesymbol, ensgid, arr[1:])
            cont = [ensgid]
            cont.append(genesymbol)
            cont.extend(arr[1:])
            f.write('\t'.join(cont) + '\n')
        except KeyError:
            try:
                ensgid = map_gene_ensgid2[genesymbol.upper()]
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

def modify_ontology_file(incsvfile, outfile):
    f = open(outfile, 'w')
    h = {}
    nomap = {}
    acronym_list  = []
    with open(incsvfile, newline='') as cf:
        for arr in csv.reader(cf, delimiter=',', quotechar='"'):
            if arr[0] == "id":
                for k in range(len(arr)):
                    h[arr[k]] = k
                arr[h['id']] = '#' + arr[h['id']]
            cont = [arr[h['id']]]

            acronym = arr[h['acronym']]
            name = arr[h['name']]
            if acronym in acronym_list:
                if 'left' in name.lower():
                    acronym = acronym + 'l'
                if 'right' in name.lower():
                    acronym = acronym + 'r'
            if acronym in acronym_list:
                try:
                    nomap[acronym] += 1
                except KeyError:
                    nomap[acronym] = 2
                acronym = acronym + str(nomap[acronym])
            acronym = arr[h['acronym']]

            # cont.append(acronym)
            # cont.append(name)

            cont = arr

            f.write('\t'.join(cont) + '\n')

            acronym_list.append(acronym)
    f.close()
    print('Saved', outfile)


def merge_expression_file(out1, out2, merged_out):
    pass


if __name__ == "__main__":
    import preproc_util
    import file_util
    path = preproc_util.DATASOURCEPATH + "/EXPRESSION/BRAIN_ATLAS/"
    out = "RNAseqTPM.tsv"
    # modify_ontology_file(path + "rnaseq_donor10021/Ontology.csv", path + "Ontology.tsv")
    preproc_BRAINATLAS(path + "rnaseq_donor10021/", "RNAseqTPM.csv", "SampleAnnot.csv", path + "Ontology.tsv", out)
    preproc_BRAINATLAS(path + "rnaseq_donor9861/", "RNAseqTPM.csv", "SampleAnnot.csv", path + "Ontology.tsv", out)
    # out1 = path + "rnaseq_donor10021/" + out
    # out2 = path + "rnaseq_donor9861/" + out
    # merged_out = path + "BrainAtlas_RNAseqTPM_10021_9861.tsv"
    # merge_expression_file(out1, out2, merged_out)
