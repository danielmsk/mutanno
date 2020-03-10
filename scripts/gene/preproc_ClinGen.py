#!/usr/bin/env python
# -*- coding: utf-8 -*-
# clingen_pars.py
# made by Daniel Minseok Kwon
# 2020-03-03 15:42:28
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
import csv

# <th style="width:22%" class="text-center">Clinical Validity Classifications</th>
#               <th style="width:22%" class="text-center">Evidence-Based Summary</th>
#               <th style="width:22%" class="text-center">Haploinsufficiency Score</th>
#               <th style="width:22%" class="text-center">Triplosensitivity Score</th>

def strip_tag(html):
    html = html.replace('</td>','').strip()
    html = str_util.strip_tag(html).strip()
    return html

def pars_clingen_gene_list(cont):
    genelist = []
    scont = str_util.substr(cont, '<tbody>', '</tbody>').strip()
    for block in scont.split('<tr>'):
        block = block.strip()
        if block != "":
            gblock = str_util.substr(block, "/kb/genes/HGNC:", "</a>")
            arr2 = gblock.split('">')
            hgnc_id = arr2[0].strip()
            gene = arr2[1].strip()


            arr2 = block.split('<td>')
            clinical_validity = strip_tag(arr2[1])
            evidence = strip_tag(arr2[2])
            haploinsufficiency_score = strip_tag(arr2[3])
            triplosensitivity_score = strip_tag(arr2[4])
            print(hgnc_id, gene)

            if hgnc_id != "":
                d = {}
                d['hgnc_id'] = hgnc_id
                d['gene'] = gene
                d['clinical_validity'] = clinical_validity
                d['evidence'] = evidence
                d['haploinsufficiency_score'] = haploinsufficiency_score
                d['triplosensitivity_score'] = triplosensitivity_score
                # genelist.append([hgnc_id, gene])
                genelist.append(d)
            # break
    return genelist

def download_list(htmlfile):
    print ("Downloading..", url)
    cont = web_util.get_url(url)
    print ("Saving..", htmlfile)
    file_util.fileSave(htmlfile, cont, 'w')
    return cont

def save_tsv(htmlfile,datalist):
    hgnc = loading_hgnc()
    out = htmlfile + ".mod.tsv"

    fieldlist = list(datalist[0].keys())

    # outcont = "#HCNC\tGENE\tensgid\tclinical_validity\tevidence\thaploinsufficiency_score\ttriplosensitivity_score\n"
    outcont = "#" + '\t'.join(fieldlist) + '\tensgid\n'
    for d in datalist:
        ensgid = ""
        for j in range(len(fieldlist)):
            if j > 0:
                outcont += '\t'
            f1 = fieldlist[j]
            outcont += d[f1]
            if f1 == "hgnc_id":
                ensgid = hgnc[d['hgnc_id']]
        outcont += '\t' + ensgid
        outcont += '\n'
    file_util.fileSave(out, outcont, 'w')
    print("Saved", out)

def loading_hgnc():
    global hgnc_file
    hgnc = {}
    h = {}
    for line in file_util.gzopen(hgnc_file):
        line = line.decode('UTF-8')
        if line[0] == "#":
            arr = line[1:].split('\t')
            for k in range(len(arr)):
                h[arr[k]] = k
        else:
            arr = line.split('\t')
            hgnc[arr[h['hgnc_id']].strip()] = arr[h['ensembl_gene_id']].strip()
    return hgnc

def preprocessing_clingen_gene(url):
    timekey = time_util.getToday()
    htmlfile = path + "clingen_list_" + timekey + ".html"

    if file_util.is_exist(htmlfile):
        cont = file_util.fileOpen(htmlfile)
    else:
        cont = download_list(htmlfile)
    
    clingen_gene_list = pars_clingen_gene_list(cont)
    print("Extracted " + str(len(clingen_gene_list)) + " genes..")

    save_tsv(htmlfile,clingen_gene_list)

def combine_disease_based_on_gene(clingen_validity):
    rstmap = {}
    for d in clingen_validity:
        try:
            for f1 in d.keys():
                rstmap[d['hgnc_id']][f1].append(d[f1])
        except KeyError:
            m = {}
            for f1 in d.keys():
                m[f1] = [d[f1]]
            rstmap[d['hgnc_id']] = m

    rstlist = []
    for hgnc_id in rstmap.keys():
        d = {}
        d['hgnc_id'] = hgnc_id
        m = rstmap[hgnc_id]
        for f1 in m.keys():
            if f1 == 'hgnc_id' or f1 == 'GENE_SYMBOL':
                d[f1] = m[f1][0]
            else:
                d[f1] = '|'.join(m[f1])
        rstlist.append(d)
    return rstlist


def preprocessing_clingen_gene_validity_curations(url):
    timekey = time_util.getToday()
    csvfile = path + "clingen_gene-validity_" + timekey + ".csv"

    if not file_util.is_exist(csvfile):
        cont = download_list(csvfile)

    no_headerline = 0
    clingen_validity = []
    with open(csvfile, newline='') as f:
        spamreader = csv.reader(f, delimiter=',', quotechar='"')
        for row in spamreader:
            if row[0] == "GENE SYMBOL":
                header = []
                for k in range(len(row)):
                    header.append(row[k].split('(')[0].strip().replace(' ','_'))
                header[header.index('GENE_ID')] = 'hgnc_id'
            elif '++++++' in row[0]:
                no_headerline += 1
            elif no_headerline >= 2:
                d = {}
                for k in range(len(row)):
                    d[header[k]] = row[k]
                d['hgnc_id'] = d['hgnc_id'].replace('HGNC:','')
                d['CLASSIFICATION_DATE'] = d['CLASSIFICATION_DATE'].replace('T',' ').replace('Z','')
                clingen_validity.append(d)
    save_tsv(csvfile,clingen_validity)

    clingen_validity_genebased = combine_disease_based_on_gene(clingen_validity)
    save_tsv(csvfile.replace('.csv','') + '_genebased.csv', clingen_validity_genebased)


if __name__ == "__main__":
    import proc_util
    import file_util
    import web_util
    import time_util
    import str_util
    url = "https://search.clinicalgenome.org/kb/curations"
    path = "/home/mk446/mutanno/DATASOURCE/GENE/CLINGEN/"
    hgnc_file = "/home/mk446/mutanno/DATASOURCE/GENE/HGNC/hgnc_complete_set.mod.txt.gz"
    preprocessing_clingen_gene(url)
    
    url = "https://search.clinicalgenome.org/kb/gene-validity.csv"
    preprocessing_clingen_gene_validity_curations(url)


