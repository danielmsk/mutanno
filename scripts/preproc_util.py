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
ENSEMBLgene = DATASOURCEPATH + "/ENSEMBL/hg38/Homo_sapiens.GRCh38.99.bed.gz"
HGNC = DATASOURCEPATH + "/GENE/HGNC/hgnc_complete_set.mod.txt.gz"


def get_map_genesymbol_ensgid_from_hgnc():
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
            m[arr[h['symbol']].strip()] = arr[h['ensembl_gene_id']]
            for asymbol in arr[h['alias_symbol']].replace('"', '').split('|'):
                m[asymbol.strip()] = arr[h['ensembl_gene_id']]
            for asymbol in arr[h['prev_symbol']].replace('"', '').split('|'):
                m[asymbol.strip()] = arr[h['ensembl_gene_id']]

    return m


def get_map_genesymbol_ensgid_from_ensembl():
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
            m[arr[h['gene_symbol']]] = arr[h['ensgid']]

    return m


if __name__ == "__main__":

    pass
