import sys
import os
import file_util
import tabix

entrezmap = {}
refseqmap = {}

def convert_uniprot_transmem(desc_value):
    arr = []
    for v1 in desc_value.split(';'):
        arr.append(v1.strip())
    return ";".join(arr)

def load_entrez_refseq_id():
    path = "/home/mk446/bio/mutanno/DATASOURCE/ENSEMBL/"
    entrezmap = {}
    for line in file_util.gzopen(path + 'Homo_sapiens.GRCh38.98.entrez.tsv.gz'):
        line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        entrezmap[arr[0].strip()] = arr[3].strip()
    refseqmap = {}
    for line in file_util.gzopen(path + 'Homo_sapiens.GRCh38.98.refseq.sorted.tsv.gz'):
        line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        enst_id = arr[1].strip()
        refseqmap[enst_id] = arr[3].strip()
    return entrezmap, refseqmap


def get_entrez_id(gene):
    global entrezmap, refseqmap
    if len(entrezmap.keys()) == 0:
        entrezmap, refseqmap = load_entrez_refseq_id()
    try:
        eid = entrezmap[gene]
    except KeyError:
        eid = ''
    return eid


def get_refseq_id(transcriptid):
    global entrezmap, refseqmap
    if len(refseqmap.keys()) == 0:
        entrezmap, refseqmap = load_entrez_refseq_id()
    try:
        tid = refseqmap[transcriptid]
    except KeyError:
        tid = ''
    return tid

