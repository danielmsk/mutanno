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
from src.mutanno.util.vcf_util import INFOIDX
from src.mutanno.util.struct_util import get_dict_value as dv
from src.mutanno.validate import AnnotVCFValidator, AnnotTSIValidator


test_region = []
# test_region.append('M:10-169')
# test_region.append('1:155197110-155197115') # CLINVAR no variant, but TAG is remained.
# test_region.append('1:59900001-59900010')
# test_region.append('1:39984-40000')
# test_region.append('4:49158516-49158518') #TOPMED
# test_region.append('19:23017580-23017580') #TOPMED
# test_region.append('2:156610760-156610775') #TOPMED
# test_region.append('2:100000-600000') #TOPMED
test_region.append('2:33814-33824') #TOPMED


def validate_annottsi(tsi, dsjson, opt):
    va = AnnotTSIValidator(tsi)
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
        # arg.extend(['-region','1:22500-22501']) # gnomAD only variant
        # arg.extend(['-region','1:952029-952609']) # clinvar variant
        # arg.extend(['-region','1:952029-952030']) # clinvar variant
        # arg.extend(['-region','1:10163-10174']) # clinvar variant
        # arg.extend(['-region','1:15034-15035']) # clinvar varian
        # arg.extend(['-region','1:10001-10100']) # clinvar variant
        # arg.extend(['-region','5:109000001-109000002']) 
        sys.argv = arg
        print_command()
        mutanno.cli()
        validate_annottsi(out, dsjson, arg)

def test_make_fullannot_sourcefile():
    ds = test_conf.get_ds()
    # dsjson = ds['fullannot']['ds']
    
    i = 0
    for r1 in test_region:
        i += 1
        key = "fullannot"
        dsjson = ds[key]['ds']
        out = "out/"+key+".datasource_test_"+str(i)+".tsi"

        arg = ['mutanno']
        arg.append('makedata')
        arg.extend(['-ds',dsjson])
        arg.extend(['-out',out])
        arg.extend(['-vartype','SNV'])
        # arg.extend(['-blocksize','5000'])
        arg.extend(['-region',r1])
        arg.extend(['-target_source','VEP'])

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
        validate_annottsi(out, dsjson, arg)


if __name__ == "__main__":
    # test_make_microannot_sourcefile()
    test_make_fullannot_sourcefile()
    



