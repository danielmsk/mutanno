import os
import sys
import time
# sys.path.append('..')
# from src.mutanno.util import file_util
# from src.mutanno.util import proc_util
import mutanno

VCF_VALIDATOR = "vcf-validator"
BGZIP = "bgzip"

TESTS_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(TESTS_DIR, 'data')
OUT_DIR = os.path.join(TESTS_DIR, 'out')
DS_DIR = os.path.join(TESTS_DIR, 'data_structure_json')
CHAINFILE = os.path.join(DATA_DIR, 'hg38ToHg19.over.chain.gz')

REPLACETAG = {}
REPLACETAG['TESTDATA_PATH'] = DATA_DIR
REPLACETAG['TESTOUT_PATH'] = OUT_DIR
REPLACETAG['MICRO_DS_FILE'] = os.path.join(DS_DIR, "datastructure_microannot_v0.4.5.json")
REPLACETAG['MICRO_SOURCE_FILE'] = os.path.join(DATA_DIR, "microannot.v0.4.4.target.mti.gz")
# REPLACETAG['MICRO_SOURCE_FILE'] = "/Users/pcaso/db/MUTANNO/DATASOURCE/MICROANNOT/v0.4.4_200614/microannot_datasource.v0.4.4_200614.tsi.gz"
REPLACETAG['FULL_DS_FILE'] = os.path.join(DS_DIR, "datastructure_fullannot_v0.4.8.json")
REPLACETAG['FULL_SOURCE_FILE'] = os.path.join(DATA_DIR, "fullannot.v0.4.8.target.mti.gz")
REPLACETAG['CHAINFILE'] = CHAINFILE

DEFAULT_VCF_VALIDATOR_MSG = "The header tag 'reference' not present. (Not required but highly recommended.)"


def getNow():
    now = time.localtime()
    s = "%04d%02d%02d_%02d%02d%02d" % (
        now.tm_year, now.tm_mon, now.tm_mday, now.tm_hour, now.tm_min, now.tm_sec)
    return s


def get_ds():
    global DS_DIR
    ds = {}
    for fname in mutanno.util.file_util.listdir(DS_DIR, '.json'):
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
        mutanno.util.file_util.save_gzip(out)
        out2 = out + '.gz'
    return mutanno.util.proc_util.run_cmd(VCF_VALIDATOR + ' ' + out2).strip()


def comp_previous_out(out, prevout):
    outcont = mutanno.util.file_util.fileOpen(out)
    prevoutcont = mutanno.util.file_util.fileOpen(prevout)
    assert outcont == prevoutcont, "Not matched: " + out + " " + prevout


def check_vcf_validator(out):
    # assert test_conf.run_vcf_validator(out) == DEFAULT_VCF_VALIDATOR_MSG
    assert run_vcf_validator(out) == ''


def validate_annotvcf(vcf, dsjson, opt):
    va = mutanno.validate.AnnotVCFValidator(vcf)
    va.set_datastructure(jsonfile=dsjson)
    va.set_option(opt)
    va.validate()


def validate_annottsi(tsi, dsjson, opt):
    va = mutanno.validate.AnnotTSIValidator(tsi)
    va.set_datastructure(jsonfile=dsjson)
    va.set_option(opt)
    va.validate()

# from src.mutanno.util.struct_util import get_dict_value as dv
# def load_dsjson(dsjson):
#     js = json.load(open(dsjson,'r'))
#     rst = {}
#     for s1 in js['source']:
#         # print(s1['name'])
#         if dv(s1, 'is_available', True):
#             for f1 in s1['fields']:
#                 fname = dv(f1, 'name2', f1['name'])
#                 k1 = s1['name'].lower() + '_' + fname.lower()
#                 f1['subembedded'] = dv(s1, 'subembedded', '')
#                 rst[k1] = f1
#     return rst
