#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import filecmp
import json
import test_conf
sys.path.append('..')
from src import mutanno
from src.mutanno.util import file_util
from src.mutanno.util import proc_util
from src.mutanno.util import seq_util
from src.mutanno.util.vcf_util import INFOIDX
from src.mutanno.util.struct_util import get_dict_value as dv
from src.mutanno.validate import AnnotVCFValidator, AnnotTSIValidator


test_region = []
# test_region.append('M:10-169')
# test_region.append('1:155197110-155197115') # CLINVAR no variant, but TAG is remained.
test_region.append('1:59900001-59900010')
# test_region.append('1:1-1000000')
# test_region.append('4:49158516-49158518') #TOPMED
# test_region.append('19:23017580-23017580') #TOPMED
# test_region.append('2:156610760-156610775') #TOPMED
# test_region.append('2:100000-600000') #TOPMED
# test_region.append('2:33814-33815') #TOPMED
# test_region.append('1:2120000-2225000') #TOPMED
# test_region.append('1:2144940-2144980')
# test_region.append('1:73120540-74181530') #TOPMED

# for chrom in seq_util.CHROM_LEN['hg38'].keys():
#     test_region.append(chrom + ':1-' + str(seq_util.CHROM_LEN['hg38'][chrom])) #chrom1 all


def validate_annottsi(tsi, dsjson, opt, sourcename=''):
    va = AnnotTSIValidator(tsi, sourcename)
    va.set_datastructure(jsonfile=dsjson)
    va.set_option(opt)
    va.validate()

def print_command():
    print('>>command: ' + ' '.join(sys.argv))

def test_make_microannot_sourcefile():
    ds = test_conf.get_ds()
    # dsjson = ds['fullannot']['ds']
    i = 0
    for r1 in test_region:
        i += 1
        key = "microannot"
        dsjson = ds[key]['dso2']
        out = "out/"+key+".datasource_test_"+str(i)+".tsi"

        arg = ['mutanno']
        arg.append('makedata')
        arg.extend(['-ds',dsjson])
        arg.extend(['-out',out])
        arg.extend(['-vartype','SNV'])
        arg.extend(['-blocksize','5000'])
        arg.extend(['-region',r1])

        sys.argv = arg
        print_command()
        mutanno.cli()
        # validate_annottsi(out, dsjson, arg)

def test_make_fullannot_sourcefile():
    ds = test_conf.get_ds()
    # dsjson = ds['fullannot']['ds']
    
    i = 0
    for r1 in test_region:
        i += 1
        key = "fullannot"
        dsjson = ds[key]['ds']
        # out = "out/"+key+".datasource_test_"+str(i)+".tsi"
        out = "/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/test_full_"+str(i)+".tsi"

        arg = ['mutanno']
        arg.append('makedata')
        arg.extend(['-ds',dsjson])
        arg.extend(['-out',out])
        arg.extend(['-vartype','SNV'])
        arg.extend(['-blocksize','10000'])
        arg.extend(['-region',r1])
        # arg.extend(['-region','all'])
        # arg.extend(['-target_source','VEP'])
        # arg.extend(['-target_source','UK10K'])
        # arg.extend(['-apply_datastructure'])

        # arg.extend(['-region','1:10001-10002']) # clinvar variant
        # arg.extend(['-region','1:22500-22501']) # gnomAD only variant
        # arg.extend(['-region','1:952029-952659']) # clinvar variant
        # arg.extend(['-region','1:952029-952030']) # clinvar variant
        # arg.extend(['-region','1:10163-20174']) # clinvar variant
        # arg.extend(['-region','1:15034-15035']) # clinvar variant
        # arg.extend(['-region','1:39984-40000']) # clinvar variant
        
        sys.argv = arg
        print_command()
        mutanno.cli()
        # validate_annottsi(out, dsjson, arg, 'UK10K')


if __name__ == "__main__":
    test_make_microannot_sourcefile()
    test_make_fullannot_sourcefile()
    



