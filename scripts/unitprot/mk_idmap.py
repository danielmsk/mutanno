#!/usr/bin/env python
# -*- coding: utf-8 -*-
# mk_idmap.py
# made by Daniel Minseok Kwon
# 2020-02-04 10:10:50
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


idmap_selected_col = []
idmap_selected_col.append("UniProtID")
idmap_selected_col.append("UniProtKB-ID")
idmap_selected_col.append("GeneID")
idmap_selected_col.append("RefSeq")
idmap_selected_col.append("GI")
idmap_selected_col.append("PDB")
idmap_selected_col.append("GO")
idmap_selected_col.append("UniRef100")
idmap_selected_col.append("UniRef90")
idmap_selected_col.append("UniRef50")
idmap_selected_col.append("UniParc")
idmap_selected_col.append("PIR")
idmap_selected_col.append("NCBI_TaxID")
idmap_selected_col.append("MIM")
idmap_selected_col.append("")
idmap_selected_col.append("PubMed1")
idmap_selected_col.append("EMBL")
idmap_selected_col.append("EMBL-CDS")
idmap_selected_col.append("Ensembl")
idmap_selected_col.append("Ensembl_TRS")
idmap_selected_col.append("Ensembl_PRO")
idmap_selected_col.append("PubMed2")

def rm_idmap_version(id):
    return id.split('-')[0]

def save_idmap(uniprotid, idmap, m):
    for k1 in m.keys():
        try:
            tmp = idmap[uniprotid]
        except KeyError:
            idmap[uniprotid] = {}
        try:
            tmp = idmap[uniprotid][k1]
        except KeyError:
            idmap[uniprotid][k1] = []

        idmap[uniprotid][k1].extend(m[k1])
        # for xid in m[k1]:
        #     if xid not in idmap[uniprotid][k1]:
        #         idmap[uniprotid][k1].append(xid)
    return idmap

def load_idmap_file(idmap):
    print("load_idmap_file")
    i = 0
    prev_uniprotid = ""
    m = {}
    for line in file_util.gzopen(idmap_file):
        i += 1
        line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        # print(arr)
        uniprotid = rm_idmap_version(arr[0].strip())
        
        if prev_uniprotid != uniprotid and len(m.keys()) > 1:
            idmap = save_idmap(uniprotid, idmap, m)
            
            m = {}

        xid = arr[2].strip()
        site = arr[1].strip()
        try:
            m[site].append(xid)
        except KeyError:
            m[site] = [xid]

        prev_uniprotid = uniprotid
        if i % 100000 == 0:
            print(i, arr)
            # break
            pass

    if len(m.keys()) > 1:
        save_idmap(uniprotid, idmap, m)

    return idmap

def load_idmap_selected_file(idmap):
    print("load_idmap_selected_file")
    i = 0
    for line in file_util.gzopen(idmap_seleted_file):
        i += 1
        line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()

        uniprotid = arr[idmap_selected_col.index('UniProtID')].strip()
        for j in range(len(idmap_selected_col)):
            site = idmap_selected_col[j]
            try:
                tmp = idmap[uniprotid]
            except KeyError:
                idmap[uniprotid] = {}
            try:
                tmp = idmap[uniprotid][site]
            except KeyError:
                idmap[uniprotid][site] = []

            idmap[uniprotid][site].extend(arr[j].strip().split('; '))
            # for xid in arr[j].strip().split('; '):
            #     xid = xid.strip()
            #     if xid not in idmap[uniprotid][site]:
            #         idmap[uniprotid][site].append(xid)

        if i % 10000 == 0:
            print(i, arr[:5])
            # break
            pass
    # print(idmap.keys())
    # print(idmap['Q00005'])
    return idmap

def load_ensembl_uniprot_map(idmap):
    print("load_ensembl_uniprot_map")
    i = 0
    for line in file_util.gzopen(ensembl_uniprot_map):
        i += 1
        line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        # print(arr)
        if arr[0] != 'gene_stable_id':
            ensgid = arr[0].strip()
            enstid = arr[1].strip()
            enspid = arr[2].strip()
            uniprotid = arr[3].strip()
            try:
                tmp = idmap[ensgid]
            except KeyError:
                idmap[ensgid] = {'Ensembl_TRS':[],'Ensembl_PRO':[],'UniProtID':[]}

            idmap[ensgid]['Ensembl_TRS'].append(enstid)
            idmap[ensgid]['Ensembl_PRO'].append(enspid)
            idmap[ensgid]['UniProtID'].append(uniprotid)

        if i % 10000 == 0:
            print(i, arr)
            # break
            pass

    # print(idmap)

    return idmap

def save_all_merged_xref_map(idmap, idmap_uniprot):
    print("save_all_merged_xref_map")

    f = open(out, 'w')
    cont = ['EnsemblID']
    cont.append('Ensembl_TRS')
    cont.append('Ensembl_PRO')
    cont.append('UniProtID')

    uniprot_site_list = []
    for uniprotid in idmap_uniprot.keys():
        for k1 in idmap_uniprot[uniprotid].keys():
            if k1 not in uniprot_site_list:
                if k1 != '':
                    uniprot_site_list.append(k1)
    cont.extend(uniprot_site_list)
    f.write('\t'.join(cont) + '\n')
    i = 0
    for ensgid in idmap.keys():
        i += 1
        cont = []
        cont.append(ensgid)
        cont.append('|'.join(idmap[ensgid]['Ensembl_TRS']))
        cont.append('|'.join(idmap[ensgid]['Ensembl_PRO']))
        cont.append('|'.join(idmap[ensgid]['UniProtID']))

        uniprotmap = {}
        for uniprotid in idmap[ensgid]['UniProtID']:
            for site in uniprot_site_list:
                try:
                    tmp = idmap_uniprot[uniprotid][site]
                    for xid in idmap_uniprot[uniprotid][site]:
                        try:
                            tmp = uniprotmap[site]
                        except KeyError:
                            uniprotmap[site] = {}
                        if xid.strip() != '' and xid != '-':
                            uniprotmap[site][xid] = 1
                except KeyError:
                    pass

        for site in uniprot_site_list:
            try:
                tmp = uniprotmap[site]
                xidlist = list(uniprotmap[site].keys())
                cont.append('|'.join(xidlist))
            except KeyError:
                cont.append('')

        if i % 1000 == 0:
            print(i)
        f.write('\t'.join(cont) + '\n')
    f.close()

    print('Saved',out)


def mk_idmap():
    
    idmap_uniprot = {}
    idmap_uniprot = load_idmap_file(idmap_uniprot)
    idmap_uniprot = load_idmap_selected_file(idmap_uniprot)
    idmap = {}
    idmap = load_ensembl_uniprot_map(idmap)
    save_all_merged_xref_map(idmap, idmap_uniprot)

    


if __name__ == "__main__":
    import proc_util
    import file_util
    out = "/home/mk446/bio/mutanno/DATASOURCE/ENSEMBL/hg38/idmap_ensembl_uniprot_xref.tsv"
    ensembl_uniprot_map = "/home/mk446/bio/mutanno/DATASOURCE/ENSEMBL/hg38/Homo_sapiens.GRCh38.98.uniprot.tsv.gz"
    idmap_file = "/home/mk446/bio/mutanno/DATASOURCE/UNIPLOT/HUMAN_9606_idmapping.dat.gz"
    idmap_seleted_file = "/home/mk446/bio/mutanno/DATASOURCE/UNIPLOT/HUMAN_9606_idmapping_selected.tab.gz"
    mk_idmap()
