#!/usr/bin/env python
# -*- coding: utf-8 -*-
# preproc_hg19_coordination.py
# made by Daniel Minseok Kwon
# 2020-04-24 01:47:05
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
sys.path.append('..')

def get_gene_id_from_gtf_format(cont):
    gene_id = ''
    for field in cont.split(';'):
        if field[:len('gene_id ')] == 'gene_id ':
            gene_id = field.replace('gene_id ','').replace('"','').strip()
            break
    return gene_id

def load_ensembl_hg19_gtf(ensembl_hg19_gtf):
    
    hg19coord = {}
    i = 0
    for line in file_util.gzopen(ensembl_hg19_gtf):
        line = file_util.decodeb(line)
        if line[0] != '#':
            i += 1
            arr = line.split('\t')

            fieldtype = arr[2]
            if fieldtype == "gene":
                chrom = arr[0]
                spos = arr[3]
                epos = arr[4]
                strand = arr[6]
                ensgid = get_gene_id_from_gtf_format(arr[8])
                if ensgid != '':
                    hg19coord[ensgid] = [chrom, spos, epos, strand]
            if i % 1000000 == 0:
                print(i)

    print ('load ' + str(len(hg19coord.keys())) + " ENSG IDs.")
    return hg19coord



def preproc_hg19_coordination(out, hg19coord):
    genesymbol2ensgid, ensgid2genesymbol = preproc_util.get_map_genesymbol_ensgid()
    # print(ensgid2genesymbol)

    f = open(out, 'w')
    cont = ['#chrom_hg19', 'spos_hg19', 'epos_hg19', 'strand_hg19','ensgid']
    f.write('\t'.join(cont) + '\n')
    ensglist = ensgid2genesymbol.keys()
    print(len(ensglist))
    cnt_miss = 0
    for ensgid in ensglist:
        try:
            coord = hg19coord[ensgid]
            cont = coord
            cont.append(ensgid)
            f.write('\t'.join(cont) + '\n')
        except KeyError:
            print(ensgid)
            cnt_miss += 1
            pass
        
    f.close()
    print('miss in hg19:', cnt_miss)
    print('Saved', out)
    

if __name__ == "__main__":
    import proc_util
    import file_util
    import preproc_util
    path = preproc_util.DATASOURCEPATH + '/'
    out = path + 'GENE/ensgid_hg19_coordination.tsv'
    ensembl_hg19_gtf = path + "ENSEMBL/hg38/Homo_sapiens.GRCh37.75.gtf.gz"

    hg19coord = load_ensembl_hg19_gtf(ensembl_hg19_gtf)
    preproc_hg19_coordination(out, hg19coord)
