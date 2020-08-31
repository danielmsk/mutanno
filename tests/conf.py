import os
import sys
import time
# sys.path.append('..')
# from src.mutanno.util import file_util
# from src.mutanno.util import proc_util
from mutanno.util import file_util
from mutanno.util import proc_util


def getNow():
    now = time.localtime()
    s = "%04d%02d%02d_%02d%02d%02d" % (
        now.tm_year, now.tm_mon, now.tm_mday, now.tm_hour, now.tm_min, now.tm_sec)
    return s


TESTS_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(TESTS_DIR, 'data')
DS_DIR = os.path.join(TESTS_DIR, 'data_structure_json')

VCF_VALIDATOR = "vcf-validator"
BGZIP = "bgzip"

REPLACETAG = {}
REPLACETAG['TESTDATA_PATH'] = DATA_DIR
REPLACETAG['MICRO_DS_FILE'] = os.path.join(DS_DIR, "datastructure_microannot_v0.4.5.json")
REPLACETAG['MICRO_SOURCE_FILE'] = os.path.join(DATA_DIR, "microannot_datasource.v0.4.4_200614_target.tsi.gz")
REPLACETAG['MICRO_SOURCE_FILE'] = "/Users/pcaso/db/MUTANNO/DATASOURCE/MICROANNOT/v0.4.4_200614/microannot_datasource.v0.4.4_200614.tsi.gz"
REPLACETAG['FULL_DS_FILE'] = os.path.join(DS_DIR, "datastructure_fullannot_v0.4.8.json")
REPLACETAG['FULL_SOURCE_FILE'] = os.path.join(DATA_DIR, "fullannot_datasource.v0.4.4_200614_target.tsi.gz")
REPLACETAG['TESTVERSION'] = getNow()

TEST_OUT1 = os.path.join(DATA_DIR, 'test_normal_snv.annot.vcf')
CHAINFILE = os.path.join(DATA_DIR, 'hg38ToHg19.over.chain.gz')

DEFAULT_VCF_VALIDATOR_MSG = "The header tag 'reference' not present. (Not required but highly recommended.)"




def get_ds():
    global DS_DIR
    ds = {}
    for fname in file_util.listdir(DS_DIR, '.json'):
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
