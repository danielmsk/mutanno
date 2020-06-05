from .util import file_util
from .util import vep_util
from .util import struct_util


def conv_consequence_delimiter(consequence):
    rst = consequence
    if isinstance(consequence, str):
        rst = consequence.replace('&', '~')
    elif isinstance(consequence, list):
        rst = []
        for c1 in consequence:
            rst.append(c1.replace('&', '~'))
    return rst


def get_score_from_vep_pathogenicity(predscore):
    rst = predscore
    if isinstance(predscore, str):
        arr = predscore.replace(')', '').split('(')
        rst = arr[1].strip()
    return rst


def get_pred_from_vep_pathogenicity(predscore):
    rst = predscore
    if isinstance(predscore, str):
        arr = predscore.replace(')', '').split('(')
        rst = arr[0].strip()
    return rst


def add_vep_most_severe(sections):
    rst = '0'
    # print (">external_functions.is_most_severe_transcript():", rst)
    most_severe = {}
    most_severe['corder'] = 9999
    target_gene = ""
    for sidx, section in enumerate(sections):
        d = {}
        # print("\tsection:",section)

        for consequence in section['Consequence']:
            if (
                (most_severe['corder'] > vep_util.VEP_CONSEQUENCE_ORDER[consequence])
                or (
                    most_severe['corder'] == vep_util.VEP_CONSEQUENCE_ORDER[consequence]
                    and section['CANONICAL'] is not None and section['CANONICAL'] == "1"
                    )
                ):
                most_severe['sidx'] = section['Gene']
                most_severe['sidx'] = sidx
                most_severe['canonical'] = section['CANONICAL']
                most_severe['corder'] = vep_util.VEP_CONSEQUENCE_ORDER[consequence]
            # elif most_severe['corder'] == vep_util.VEP_CONSEQUENCE_ORDER[consequence]:
            #     if section['CANONICAL'] is not None and section['CANONICAL'] == "1":
            # print("\t\tconsequence:",consequence)
        # print("\tmost_severe:",most_severe)

    rst = []
    ms = ''
    for sidx, section in enumerate(sections):
        if sidx == most_severe['sidx']:
            # print('====>', sections[sidx])
            sections[sidx]['MOST_SEVERE'] = '1'
            # ms = '1'
        else:
            sections[sidx]['MOST_SEVERE'] = '0'
            # ms = '0'
    # print(sections)
    # return ms
    return sections

def vep_select_biotype_add_most_severe(sections, select_biotype):
    sections = vep_select_biotype(sections, select_biotype)
    sections = add_vep_most_severe(sections)
    return sections
    

def vep_select_biotype(sections, select_biotype):
    # print (">external_functions.vep_select_biotype():", select_biotype)
    # print(sections)
    selected = []
    for sidx, section in enumerate(sections):
        # print("section['BIOTYPE']:",section['BIOTYPE'])
        if section['BIOTYPE'] in select_biotype:
            selected.append(section)
    return selected

def trim_DIP_ID(v1):
    # ex) DIP-39616N;
    rst = v1
    if isinstance(v1, str):
        rst = v1.split('-')[-1].replace('N', '')
    return rst


def remove_pdb_subversion(id_list):
    rst = id_list
    if isinstance(id_list, str):
        pidmap = {}
        for pid in id_list.split('|'):
            pidmap[pid.split(':')[0]] = 1
        pidlist = list(pidmap.keys())
        rst='|'.join(pidlist)
    return 


def convert_rmsk_strand(v1):
    rst = '0'
    if isinstance(v1, str):
        if v1 == '+':
            rst = '1'
    return rst


def convert_NS2blank(v1):
    rst = v1
    if isinstance(v1, str):
        if v1 == 'NS':
            rst = ''
    return rst


def convert_vep_strand(v1):
    rst = v1
    if isinstance(v1, str):
        rst = v1.replace('-1', '0')
    return rst


def convert_high_inf_pos(v1):
    rst = v1
    if isinstance(v1, str):
        v1.replace('Y', '1').replace('N', '0')
    return rst


def convert_canonical2boolean(canonical):
    rst = '0'
    if isinstance(canonical, str):
        if canonical == 'YES':
            rst = '1'
    return rst


def convert_uniprot_transmem(desc_value):
    rst = desc_value
    if isinstance(desc_value, str):
        arr = []
        for v1 in desc_value.split(';'):
            arr.append(v1.strip())
        rst = ";".join(arr)
    return rst


def add_genes_severe_consequence(annotdata, vcfinfo):
    # print(">external_function.add_severe_consequence()")
    # print('\t',annotdata['VEP'])

    vep_sections = []
    if annotdata is not None and 'VEP' in annotdata.keys():
        vep_sections = annotdata['VEP']
    elif vcfinfo is not None and 'VEP' in vcfinfo.keys():
        vep_sections = vcfinfo['VEP']

    for attr in vep_sections:
        if attr['MOST_SEVERE'] == '1':
            ensg = attr['Gene']
            most_severe_transcript = attr['Feature']
            most_severe_consequence = attr['Consequence']
            break
    if len(vep_sections) == 0:
        rst_sections = []
    else:
        rst_sections =[{'ensg':ensg, 'most_severe_transcript':most_severe_transcript, 'most_severe_consequence':most_severe_consequence}]
    return rst_sections


# entrezmap = {}
# refseqmap = {}
# def load_entrez_refseq_id():
#     path = "/home/mk446/bio/mutanno/DATASOURCE/ENSEMBL/hg38/"
#     entrezmap = {}
#     for line in file_util.gzopen(path + 'Homo_sapiens.GRCh38.98.entrez.tsv.gz'):
#         line = line.decode('UTF-8')
#         arr = line.split('\t')
#         arr[-1] = arr[-1].strip()
#         entrezmap[arr[0].strip()] = arr[3].strip()
#     refseqmap = {}
#     for line in file_util.gzopen(path + 'Homo_sapiens.GRCh38.98.refseq.sorted.tsv.gz'):
#         line = line.decode('UTF-8')
#         arr = line.split('\t')
#         arr[-1] = arr[-1].strip()
#         enst_id = arr[1].strip()
#         refseqmap[enst_id] = arr[3].strip()
#     return entrezmap, refseqmap


# def get_entrez_id(gene):
#     global entrezmap, refseqmap
#     if len(entrezmap.keys()) == 0:
#         entrezmap, refseqmap = load_entrez_refseq_id()
#     try:
#         eid = entrezmap[gene]
#     except KeyError:
#         eid = ''
#     return eid


# def get_refseq_id(transcriptid):
#     global entrezmap, refseqmap
#     if len(refseqmap.keys()) == 0:
#         entrezmap, refseqmap = load_entrez_refseq_id()
#     try:
#         tid = refseqmap[transcriptid]
#     except KeyError:
#         tid = ''
#     return tid
