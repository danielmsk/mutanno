import sys
import shlex
import filecmp
import json
import time
import conf

import mutanno
# sys.path.append('..')
# from src import mutanno

prog = "mutanno"

MICROANNOT_CMD = """
    annot \
    -vcf #VCF# \
    -ds  #MICRO_DS_FILE# \
    -out #OUT# \
    -sourcefile #MICRO_SOURCE_FILE# \
    -split_multi_allelic_variant \
    -genoinfo \
    -use_raw_source
"""

FULLANNOT_CMD = """
    annot \
    -vcf #VCF# \
    -ds  #FULL_DS_FILE# \
    -out #OUT# \
    -sourcefile #FULL_SOURCE_FILE# \
    -hg19 \
    -chain #CHAINFILE# \
    -clean_tag MUTANNO SpliceAI CLINVAR gnomADgenome
"""

TEST_VERSION = conf.getNow()

def test_microannot():
    cmd = MICROANNOT_CMD
    opt={}
    opt['VCF'] = '#TESTDATA_PATH#/TRIO_GAPFIS8ZSPEO.target.vcf.gz'
    opt['PREV_OUT'] = '#TESTDATA_PATH#/TRIO_GAPFIS8ZSPEO.target.vcf.gz.microannot_v0.4.4.vcf.gz'
    opt['OUT'] = '#TESTOUT_PATH#/TRIO_GAPFIS8ZSPEO.target.vcf.gz.microannot_v0.4.4.'+TEST_VERSION+'.vcf'
    optkeys = list(opt.keys())
    for t1 in optkeys:
        for t2 in conf.REPLACETAG.keys():
            if '#'+t2+'#' in opt[t1]:
                opt[t1] = opt[t1].replace('#'+t2+'#', conf.REPLACETAG[t2])
            if t2 not in opt.keys():
                opt[t2] = conf.REPLACETAG[t2]

    for t1 in opt.keys():
        cmd = cmd.replace('#'+t1+'#', opt[t1])
    cmd = prog + " " + cmd.strip()
    sys.argv = shlex.split(cmd)
    print(' '.join(sys.argv))
    # print(shlex.quote(sys.argv))
    mutanno.cli()
    conf.check_vcf_validator(opt['OUT'])
    conf.comp_previous_out(opt['OUT'], opt['PREV_OUT'])
    conf.validate_annotvcf(opt['OUT'], opt['MICRO_DS_FILE'], sys.argv)

def test_fullannot():
    cmd = FULLANNOT_CMD
    opt = {}
    opt['VCF'] = '#TESTOUT_PATH#/TRIO_GAPFIS8ZSPEO.target.vcf.gz.microannot_v0.4.4.'+TEST_VERSION+'.vcf'
    opt['PREV_OUT'] = '#TESTDATA_PATH#/TRIO_GAPFIS8ZSPEO.target.vcf.gz.microannot_v0.4.4.vcf.gz'
    opt['OUT'] = '#TESTOUT_PATH#/TRIO_GAPFIS8ZSPEO.target.vcf.gz.microannot_v0.4.4.fullannot_v0.4.8.'+TEST_VERSION+'.vcf'
    optkeys = list(opt.keys())
    for t1 in optkeys:
        for t2 in conf.REPLACETAG.keys():
            if '#'+t2+'#' in opt[t1]:
                opt[t1] = opt[t1].replace('#'+t2+'#', conf.REPLACETAG[t2])
            if t2 not in opt.keys():
                opt[t2] = conf.REPLACETAG[t2]

    for t1 in opt.keys():
        cmd = cmd.replace('#'+t1+'#', opt[t1])
    cmd = prog + " " + cmd.strip()
    sys.argv = shlex.split(cmd)
    print(' '.join(sys.argv))
    # print(shlex.quote(sys.argv))
    mutanno.cli()
    conf.check_vcf_validator(opt['OUT'])
    # conf.comp_previous_out(opt['OUT'], opt['PREV_OUT'])
    # conf.validate_annotvcf(opt['OUT'], opt['MICRO_DS_FILE'], sys.argv)


if __name__ == "__main__":
    test_microannot()
    test_fullannot()
    



# 
