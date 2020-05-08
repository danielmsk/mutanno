import file_util
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



CURRENTPATH = os.path.abspath(os.path.dirname(__file__))
DATASOURCEPATH = os.path.abspath(os.path.join(CURRENTPATH, "../../DATASOURCE/"))
TESTPATH = os.path.abspath(os.path.join(CURRENTPATH, "../tests/"))
ENSEMBLgene = DATASOURCEPATH + "/ENSEMBL/hg38/Homo_sapiens.GRCh38.99.bed.gz"
HGNC = DATASOURCEPATH + "/GENE/HGNC/hgnc_complete_set.mod.txt.gz"
NCBIgene = DATASOURCEPATH + "/NCBI/GENE/Homo_sapiens.gene_info.gz"
GENE_SYMBOL_MAP = DATASOURCEPATH + "/GENE/gene_symbolmap.tsv"

def get_symbol(symbol, flag_upper=True):
    symbol = symbol.strip()
    if flag_upper:
        symbol = symbol.upper()
    return symbol

def get_map_genesymbol_ensgid(flag_upper=True):
    global GENE_SYMBOL_MAP
    genesymbol2ensgid = {}
    ensgid2genesymbol = {}
    for line in file_util.gzopen(GENE_SYMBOL_MAP):
        line = file_util.decodeb(line)
        arr = line.split('\t')
        arr[-1] = arr[-1]
        ensgid = arr[0].strip()

        genesymbol2ensgid[get_symbol(arr[1], flag_upper)] = ensgid
        ensgid2genesymbol[ensgid] = [get_symbol(arr[1], flag_upper)]

        for s1 in arr[2].split('|'):
            if s1 not in ensgid2genesymbol[ensgid]:
                ensgid2genesymbol[ensgid].append(get_symbol(s1, flag_upper))
            genesymbol2ensgid[get_symbol(s1, flag_upper)] = ensgid

    return genesymbol2ensgid, ensgid2genesymbol


def get_map_genesymbol_ensgid_from_hgnc(flag_upper=True):
    global HGNC
    m = {}
    h = {}
    for line in file_util.gzopen(HGNC):
        line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1]
        if line[0] == "#":
            arr[0] = arr[0][1:]
            for k in range(len(arr)):
                h[arr[k]] = k
        else:
            m[get_symbol(arr[h['symbol']], flag_upper)] = arr[h['ensembl_gene_id']]
            for asymbol in arr[h['alias_symbol']].replace('"', '').split('|'):
                m[get_symbol(asymbol, flag_upper)] = arr[h['ensembl_gene_id']]
            for asymbol in arr[h['prev_symbol']].replace('"', '').split('|'):
                m[get_symbol(asymbol, flag_upper)] = arr[h['ensembl_gene_id']]

    return m


def get_map_genesymbol_ensgid_from_ensembl(flag_upper = True):
    global ENSEMBLgene
    m = {}
    h = {}
    for line in file_util.gzopen(ENSEMBLgene):
        line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1]
        if line[0] == "#":
            arr[0] = arr[0][1:]
            for k in range(len(arr)):
                h[arr[k]] = k
        else:
            m[get_symbol(arr[h['gene_symbol']], flag_upper)] = arr[h['ensgid']]

    return m


if __name__ == "__main__":

    pass
