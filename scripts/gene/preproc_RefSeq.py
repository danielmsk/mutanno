#!/usr/bin/env python
# -*- coding: utf-8 -*-
# preproc_RefSeq.py
# made by Daniel Minseok Kwon
# 2020-03-08 21:59:44
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

class ParsGbff():
    def __init__(self):
        pass
        
    def pars_gbff(self, gbff):
        print(gbff)
        self.file = gbff
        self.f = file_util.gzopen(self.file)

    def get_field(self):
        i = 0
        flag = True
        k2 = ""
        m = {}
        m['transcript_id']=[]
        m['GeneID']=[]
        m['HGNC']=[]
        noblankline = 0
        while True:
            i += 1
            line = self.f.readline()
            line = line.decode('UTF-8')

            if line.strip() == "//":
                break
            if line.strip() == "":
                noblankline += 1
            else:
                noblankline = 0
            if noblankline > 5:
                break

            k1 = line[:12].strip()
            v1 = line[12:].strip()
            # v1[-1] = v1[-1].strip()

            if k1 != "":
                k2 = k1

            if k2 == "FEATURES":
                flag = False
            if k2 == "mRNA":

                if "/transcript_id=" in v1:
                    tid = v1.replace("/transcript_id=","").replace('"','')
                    m['transcript_id'].append(tid)
                if '/db_xref="GeneID:' in v1:
                    tid = v1.replace('/db_xref="GeneID:','').replace('"','')
                    m['GeneID'].append(tid)
                if '/db_xref="HGNC:' in v1:
                    tid = v1.replace('/db_xref="HGNC:HGNC:','').replace('"','')
                    m['HGNC'].append(tid)
                # if '/db_xref="MIM:' in v1:
                #     tid = v1.replace('/db_xref="MIM:','').replace('"','')
                #     try:
                #         m['MIM'].append(tid)
                #     except KeyError:
                #         m['MIM'] = [tid]
                

            if flag:
                try:
                    m[k2] += " " + v1
                except KeyError:
                    m[k2] = v1

                # print(k1, v1)
                if i % 10000 == 0:
                    print(i, v1)

        m = self.pars_comment(m)
        return m

    def pars_comment(self, field):
        if "COMMENT" in field.keys():
            arr = field['COMMENT'].split('  Summary:')
            field['REVIEWED_REFSEQ'] = arr[0].replace('REVIEWED REFSEQ:','').strip()
            field['SUMMARY'] = ""
            if len(arr) > 1:
                field['SUMMARY'] = arr[1].strip()
        return field

def load_hgnc_ensg():
    global hgnc_file

    m = {}
    h = {}
    for line in file_util.gzopen(hgnc_file):
        line = line.decode('UTF-8')
        arr = line.split('\t')
        if line[0] == "#":
            arr[0] = arr[0][1:]
            for k in range(len(arr)):
                h[arr[k]] = k
        else:
            m[arr[h['hgnc_id']]] = arr[h['ensembl_gene_id']]
    return m

def preproc_RefSeq(gbff_file, out):
    hgnc_ensg_map = load_hgnc_ensg()

    f = open(out, 'w')
    cont = []
    cont.append('#AC')
    cont.append('engs_id')
    cont.append('transcript_id')
    cont.append('GeneID')
    cont.append('HGNC')
    # cont.append('MIM')
    cont.append('SUMMARY')
    cont.append('REVIEWED_REFSEQ')
    f.write('\t'.join(cont) + '\n')
    pg = ParsGbff()
    ensg_list = []
    for k in range(1,6):
        pg.pars_gbff(gbff_file.replace('#NO#', str(k)))
        # pg.pars_gbff(path + "aa.gz")
        i = 0
        while True:
            i += 1
            field = pg.get_field()
            if not 'ACCESSION' in field.keys():
                break

            print(field['ACCESSION'])
            # print(field)

            for j in range(len(field['transcript_id'])):

                try:
                    ensg_id = hgnc_ensg_map[field['HGNC'][j]]
                except IndexError:
                    ensg_id = ''
                if not ensg_id in ensg_list and ensg_id != '':
                    ensg_list.append(ensg_id)
                    cont = []
                    cont.append(field['ACCESSION'])
                    cont.append(ensg_id)
                    try:
                        cont.append(field['transcript_id'][j])
                    except IndexError:
                        cont.append("")
                    try:
                        cont.append(field['GeneID'][j])
                    except IndexError:
                        cont.append("")
                    try:
                        cont.append(field['HGNC'][j])
                    except IndexError:
                        cont.append("")
                    # cont.append(field['MIM'][j])
                    cont.append(field['SUMMARY'])
                    cont.append(field['REVIEWED_REFSEQ'])
                    # print(cont)
                    f.write('\t'.join(cont) + '\n')
                
        
        # break
    f.close()
    print('Saved', out)
        

if __name__ == "__main__":
    import proc_util
    import file_util
    path = "/home/mk446/mutanno/DATASOURCE/NCBI/REFSEQ/"
    gbff_file = path + "refseqgene.#NO#.genomic.gbff.gz"
    out = path + "refseqgene.genomic.tsv"

    hgnc_file = "/home/mk446/mutanno/DATASOURCE/GENE/HGNC/hgnc_complete_set.mod.txt.gz"
    # "/home/mk446/mutanno/DATASOURCE/ENSEMBL/hg38/idmap_ensembl_uniprot_xref.tsv.gz"
    # "/home/mk446/mutanno/DATASOURCE/ENSEMBL/hg38/Homo_sapiens.GRCh38.98.entrez.tsv.gz"
    preproc_RefSeq(gbff_file, out)
