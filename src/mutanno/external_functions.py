from .util import file_util
from .util import vep_util

entrezmap = {}
refseqmap = {}


def is_most_severe_transcript(sections, section_idx, colheader):
    rst = '0'
    most_severe = {}
    target_gene = ""
    for sidx, section in enumerate(sections):
        d = {}
        d['gene'] = section[colheader.index('Gene')]
        d['canonical'] = section[colheader.index('CANONICAL')]
        d['sidx'] = sidx
        d['corder'] = 9999  ## consequence_order
        for consequence in section[colheader.index('Consequence')].split('&'):
            if d['corder'] < vep_util.VEP_CONSEQUENCE_ORDER[consequence]:
                d['corder'] = vep_util.VEP_CONSEQUENCE_ORDER[consequence]
        try:
            if most_severe[d['gene']]['corder'] > d['corder']:
                most_severe[d['gene']] = d
            elif most_severe[d['gene']]['corder'] == d['corder']:
                if most_severe[d['gene']]['canonical'] != "YES" and d['canonical'] == "YES":
                    most_severe[d['gene']] = d
        except KeyError:
            most_severe[d['gene']] = d

        if sidx == section_idx:
            target_gene = d['gene']

    if most_severe[target_gene]['sidx'] == section_idx:
        rst = '1'

    return rst



def trim_DIP_ID(v1):
    # DIP-39616N;
    v1 = v1.split('-')[-1].replace('N', '')
    return v1


def remove_pdb_subversion(id_list):
    pidmap = {}
    for pid in id_list.split('|'):
        pidmap[pid.split(':')[0]] = 1
    pidlist = list(pidmap.keys())
    return '|'.join(pidlist)


def convert_rmsk_strand(v1):
    r1 = '0'
    if v1 == '+':
        r1 = '1'
    return r1


def convert_NS2blank(v1):
    if v1 == 'NS':
        v1 = ''
    return v1


def convert_vep_strand(v1):
    return v1.replace('-1', '0')


def convert_high_inf_pos(v1):
    return v1.replace('Y', '1').replace('N', '0')


def convert_canonical2boolean(canonical):
    r1 = '0'
    if canonical == 'YES':
        r1 = '1'
    return r1


def convert_uniprot_transmem(desc_value):
    arr = []
    for v1 in desc_value.split(';'):
        arr.append(v1.strip())
    return ";".join(arr)


def load_entrez_refseq_id():
    path = "/home/mk446/bio/mutanno/DATASOURCE/ENSEMBL/hg38/"
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
