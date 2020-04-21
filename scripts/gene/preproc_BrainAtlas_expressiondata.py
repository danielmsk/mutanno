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

def load_ontology(csvfile):
    m = {}
    h = {}
    with open(csvfile, newline='') as cf:
        for arr in csv.reader(cf, delimiter=',', quotechar='"'):
            if arr[0] == "id":
                for k in range(len(arr)):
                    h[arr[k]] = k
            else:
                # m[arr[h['id']]] = arr[h['name']]
                m[arr[h['id']]] = arr[h['id']] + "." + arr[h['name']]
                
    return m

def load_ontology_tsv(tsvfile):
    m = {}
    for line in open(tsvfile):
        if line[0] != '#':
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            # m[arr[0]] = arr[1]
            m[arr[0]] = arr[0] + "." + arr[1]

    return m

def load_sampleannot(csvfile, ontology, col1="structure_id", tcol="structure_id"):
    m = []
    h = {}
    acronym_list = []
    no_map = {}
    for line in file_util.gzopen(csvfile):
        line = file_util.decodeb(line)
        arr = line.split(',')
        arr[-1] = arr[-1].strip()
        # print(arr)
        if arr[0] == col1:
            for k in range(len(arr)):
                h[arr[k]] = k
        else:
            acronym = ontology[arr[h[tcol]]]
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

def load_probe_tsv(tsvfile, genesymbol2ensgid):
    m = {}
    h = {}
    for line in file_util.gzopen(tsvfile):
        line = file_util.decodeb(line)
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        if line[0] == "#":
            arr[0] = arr[0][1:]
            for k in range(len(arr)):
                h[arr[k]] = k
        else:
            ensgid = ""
            try:
                ensgid = genesymbol2ensgid[arr[h['gene_symbol']].upper()]
                # print(ensgid)
            except KeyError:
                # print ('KeyError:', arr)
                pass
            pid = arr[h['probe_id']]
            m[pid] = {}
            m[pid]['probe_name'] = arr[h['probe_name']]
            m[pid]['gene_id'] = arr[h['gene_id']]
            m[pid]['gene_symbol'] = arr[h['gene_symbol']]
            m[pid]['ensgid'] = ensgid
        
    return m



def preproc_BRAINATLAS(dpath, tpm_csv, sampleannot_csv, probe_tsv,ontology_tsv, out):

    genesymbol2ensgid, ensgid2genesymbol = preproc_util.get_map_genesymbol_ensgid()
    

    f = open(dpath + out, 'w')
    ontology = load_ontology_tsv(ontology_tsv)
    sampleannot = load_sampleannot(dpath + sampleannot_csv, ontology)
    probemap = load_probe_tsv(dpath + probe_tsv, genesymbol2ensgid)

    cont = []
    cont.append('#probe_id')
    cont.append('probe_name')
    cont.append('gene_id')
    cont.append('gene_symbol')
    cont.append('ensgid')
    cont.extend(sampleannot)
    f.write('\t'.join(cont) + '\n')

    
    for line in file_util.gzopen(dpath + tpm_csv):
        line = file_util.decodeb(line)
        arr = line.split(',')
        arr[-1] = arr[-1].strip()
        probe_id = arr[0].replace('"', '')
        # print(arr[:5])
        try:
            probeinfo = probemap[probe_id]
            if probeinfo['ensgid'] != '':
                cont = [probe_id, probeinfo['probe_name'], probeinfo['gene_id'], probeinfo['gene_symbol'],probeinfo['ensgid']]
                cont.extend(arr[1:])
                f.write('\t'.join(cont) + '\n')
                # print(probeinfo)
        except KeyError:
            pass
            # print(genesymbol)
        # break
    # print(map_gene_ensid)
    f.close()
    print('Saved', dpath + out)



def preproc_BRAINATLAS_rnaseq(dpath, tpm_csv, sampleannot_csv, ontology_tsv, out):

    genesymbol2ensgid, ensgid2genesymbol = preproc_util.get_map_genesymbol_ensgid()

    f = open(dpath + out, 'w')
    ontology = load_ontology_tsv(ontology_tsv)
    sampleannot = load_sampleannot(dpath + sampleannot_csv, ontology, "RNAseq_sample_name","ontology_structure_id")

    cont = []
    cont.append('#ensgid')
    cont.append('gene_symbol')
    cont.extend(sampleannot)
    f.write('\t'.join(cont) + '\n')

    
    for line in file_util.gzopen(dpath + tpm_csv):
        line = file_util.decodeb(line)
        arr = line.split(',')
        arr[-1] = arr[-1].strip()
        genesymbol = arr[0].replace('"', '')
        try:
            ensgid = genesymbol2ensgid[genesymbol.upper()]
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

def get_sample_ids(expfile, sidx):
    sample_ids = []
    for line in file_util.gzopen(expfile):
        line = file_util.decodeb(line)
        arr = line[1:].split('\t')
        arr[-1] = arr[-1].strip()
        sample_ids = arr[sidx:]
        break
    return sample_ids

def merge_expression_file(samplelist, merged_out):
    expfiles = {}
    total_sample_id_list = []
    exp = {}
    for sid in samplelist:
        if 'microarray' in sid:
            out = "MicroarrayExpression.tsv.gz"
            sidx = 4
        else:
            out = "RNAseqTPM.tsv.gz"
            sidx = 2
        expfile = path + sid+"/" + out
        expfiles[sid] = expfile
        # sample_ids = get_sample_ids(expfile, sidx)

        print(expfile)
        sample_ids = []
        i = 0
        h = {}
        for line in file_util.gzopen(expfile):
            line = file_util.decodeb(line)
            arr = line[1:].split('\t')
            arr[-1] = arr[-1].strip()

            if i == 0:
                header = arr
                sample_ids = arr[sidx:]
                for k in range(len(arr)):
                    h[arr[k]] = k
            else:
                ensgid = arr[h['ensgid']]
                try:
                    exp[ensgid]
                except KeyError:
                    exp[ensgid] = {}
                # print(arr[:5])
                try:
                    exp[ensgid][sid]
                except KeyError:
                    exp[ensgid][sid] = []


                d = {}
                try:
                    d['probe_id'] = arr[h['probe_id']]
                except KeyError:
                    d['probe_id'] = ""

                for k in range(sidx, len(arr)):
                    try:
                        # if some samples have multiple expression data.
                        d[header[k]] += '~' + arr[k]
                    except KeyError:
                        d[header[k]] = arr[k]
                exp[ensgid][sid].append(d)

            i += 1

            if i % 5000 == 0:
                print(i)
                # break


        for sample_id in sample_ids:
            if sample_id not in total_sample_id_list and sample_id != 'ensgid':
                total_sample_id_list.append(sample_id)

    total_sample_id_list = sorted(total_sample_id_list)

    f = open(merged_out, 'w')

    cont = ['#ensgid']
    h1 = 'sample_id|experiment_method|probe_id'
    for sample_id in total_sample_id_list:
        h1 += '|' + sample_id.replace(' ', '_')
    cont.append(h1)
    f.write('\t'.join(cont) + '\n')

    for exsgid in exp.keys():
        cont = [exsgid]
        cont2 = []
        for sid in exp[exsgid].keys():
            sidarr = sid.split('.')

            for d in exp[exsgid][sid]:
                cont3 = []
                cont3.append(sidarr[0] + '.' + sidarr[1])
                cont3.append(sidarr[2])
                cont3.append(d['probe_id'])
                for sample_id in total_sample_id_list:
                    try:
                        cont3.append(d[sample_id])
                    except KeyError:
                        cont3.append('')

            cont2.append('|'.join(cont3))
        cont.append(','.join(cont2))
        f.write('\t'.join(cont) + '\n')
    f.close()



    print(len(total_sample_id_list))
    print('Saved', merged_out)


def load_ontology_tsv_allinfo(tsvfile):
    m = {}
    for line in open(tsvfile):
        if line[0] != '#':
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            m[arr[0]] = arr

    return m

def save_structure_file(merged_out, ontology_tsv):
    ontology = load_ontology_tsv_allinfo(ontology_tsv)

    out = merged_out + '.struct.txt'
    f = open(out, 'w')
    for line in file_util.gzopen(merged_out + '.gz'):
        line = file_util.decodeb(line)
        if line[0] == '#':
            ensgid = line.split('\t')[0]
            arr = line.split('\t')[1].strip().split('|')

            print(arr)
            for sid in arr[3:]:
                sidx = sid.split('.')[0]
                # cont = [sid]
                # cont.extend(ontology[sidx])
                ot = ontology[sidx]
                cont = [ot[0]]
                cont.append(ot[2] + ' (' + ot[1] + ')')
                cont.append(ot[0] + ' ' + ot[2] + ' (' + ot[1] + ') ['+ot[7]+']')
                # cont = ['{"name":"'+sid+'","title":"'+ontology[sidx][2]+'","desc":"'+ontology[sidx][0]+' '+ontology[sidx][2]+'('+ontology[sidx][1]+') [' + ontology[sidx][7]+']","is_list":true,"delimiter":"~"},']
                f.write('\t'.join(cont) + '\n')
        break

    f.close()
    print('Saved', out)

    

if __name__ == "__main__":
    import preproc_util
    import file_util
    path = preproc_util.DATASOURCEPATH + "/EXPRESSION/BRAIN_ATLAS/"
    ontology_tsv = path + "Ontology.tsv"
    
    # modify_ontology_file(path + "rnaseq_donor10021/Ontology.csv", path + "Ontology.tsv")
    samplelist = []
    samplelist.append("H0351.1009.microarray")
    samplelist.append("H0351.1012.microarray")
    samplelist.append("H0351.1015.microarray")
    samplelist.append("H0351.1016.microarray")
    samplelist.append("H0351.2001.microarray")
    samplelist.append("H0351.2002.microarray")
    samplelist.append("H0351.2001.rnaseq")
    samplelist.append("H0351.2002.rnaseq")



    for sid in samplelist:
        if 'microarray' in sid:
            out = "MicroarrayExpression.tsv"
            # preproc_BRAINATLAS(path + ""+sid+"/", "MicroarrayExpression.csv.gz", "SampleAnnot.csv.gz", "Probes.tsv.gz" ,path + "Ontology.tsv", out)
        else:
            out = "RNAseqTPM.tsv"
            # preproc_BRAINATLAS_rnaseq(path + sid+"/", "RNAseqTPM.csv.gz", "SampleAnnot.csv.gz", path + "Ontology.tsv", out)
            

    merged_out = path + "BrainAtlas_rnaseq_microarray.tsv"
    merge_expression_file(samplelist, merged_out)
    # save_structure_file(merged_out, ontology_tsv)

    
