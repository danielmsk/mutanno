import sys
import filecmp
import json
import test_conf
sys.path.append('..')
from src import mutanno
from src.mutanno.util import file_util
from src.mutanno.util import proc_util
from src.mutanno.util.vcf_util import INFOIDX
from src.mutanno.util.struct_util import get_dict_value as dv
from src.mutanno.validate import AnnotVCFValidator, AnnotTSIValidator


def test_annotvcf_veiwer():
    ds = test_conf.get_ds()

    # vcf = "/home/mk446/mutanno/TEST/0616/GAPFIP83PL7E_test2.vcf.gz.fullannot.vcf.gz"
    # vcf = "/home/mk446/mutanno/TEST/0616/GAPFI3JX5D2J_test1.vcf.gz"
    # vcf = "/home/mk446/mutanno/TEST/0616/GAPFI3JX5D2J_test1.vcf.gz.fullannot.vcf.gz"
    vcf = "/home/mk446/mutanno/TEST/0616/GAPFI3JX5D2J_test2.vcf.gz.fullannot.vcf.gz"
    dsjson = ds['fullannot']['clean']

    arg = ['mutanno']
    arg.append('view')
    arg.extend(['-vcf',vcf])
    arg.extend(['-ds',dsjson])
    # arg.extend(['-region',"chr1:30663-30663"])
    # arg.extend(['-region',"chr22:37971233-37971233"])
    arg.extend(['-region',"chr1:1041200-1041200"])

    # ENST00000379370
    
    sys.argv = arg
            
    print('>>command:',' '.join(sys.argv))
    mutanno.cli()
    # check_vcf_validator(out)
    # validate_annotvcf(out, dsjson,arg)

if __name__ == "__main__":
    test_annotvcf_veiwer()

