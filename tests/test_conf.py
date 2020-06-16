import os
import sys
sys.path.append('..')

from src.mutanno.util import file_util
from src.mutanno.util import proc_util

TESTS_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(TESTS_DIR, 'data')

VCF_VALIDATOR = "/n/app/vcftools/0.1.16/bin/vcf-validator"
BGZIP = "/home/mk446/anaconda3/bin/bgzip"

SOURCEFILE = {}
# SOURCEFILE['microannot'] = os.path.join(DATA_DIR, 'mvp_datasource_v0.3.2_200309.chr#CHROM#.test.tsi.gz')
# SOURCEFILE['microannot'] = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/v0.4.2_200512/microannot_datasource.v0.4.2_20200512.tsi.gz"
# SOURCEFILE['microannot'] = "/home/mk446/bio/mutanno/DATASOURCE/MICROANNOT/v0.4.3_200604_chrom/microannot_datasource.#CHROM#.v0.4.3_200604.tsi.gz"
# SOURCEFILE['microannot'] = "/home/mk446/bio/mutanno/DATASOURCE/MICROANNOT/v0.4.3_200604/GAPFI455U3JI.tsi.gz"
SOURCEFILE['microannot'] = "/home/mk446/bio/mutanno/DATASOURCE/MICROANNOT/microannot_datasource.v0.4.4_200614.tsi.gz"

# SOURCEFILE['fullannot'] = os.path.join(DATA_DIR, 'mvp_datasource_v0.3.2_200309.chr#CHROM#.test.tsi.gz')
# SOURCEFILE['fullannot'] = os.path.join(DATA_DIR, 'mvp_datasource_v0.3.2_200309.chr#CHROM#.test.tsi.gz')
SOURCEFILE['fullannot'] = '/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/v0.4.0_200608_chrom/mvp_datasource_v0.4.0_200608.chr#CHROM#.tsi.gz'
# SOURCEFILE['GENE'] 

TEST_OUT1 = os.path.join(DATA_DIR, 'test_normal_snv.annot.vcf')
CHAINFILE = os.path.join(DATA_DIR, 'hg38ToHg19.over.chain.gz')

DEFAULT_VCF_VALIDATOR_MSG = "The header tag 'reference' not present. (Not required but highly recommended.)"


def get_ds():
    global DATA_DIR
    ds = {}
    for fname in file_util.listdir(DATA_DIR, '.json'):
        if fname.startswith('datastructure_'):
            arr = fname.split('_')

            try:
                ds[arr[1]]
            except KeyError:
                ds[arr[1]] = {}
            
            if 'dso2' in arr[2]:
                ds[arr[1]]['dso2'] = os.path.join(DATA_DIR, fname)
            elif 'ds' in arr[2]:
                ds[arr[1]]['ds'] = os.path.join(DATA_DIR, fname)
            else:
                ds[arr[1]]['clean'] = os.path.join(DATA_DIR, fname)
    
    return ds

def run_vcf_validator(out):
    global VCF_VALIDATOR
    out2 = out
    if not out.endswith('.gz'):
        file_util.save_gzip(out)
        out2 = out + '.gz'
    return proc_util.run_cmd(VCF_VALIDATOR + ' ' + out2).strip()