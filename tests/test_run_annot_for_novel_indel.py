import sys
import shlex
import filecmp
import json
import time
import conf

# import mutanno
sys.path.append('..')
from src import mutanno

prog = "mutanno"

PREPROCESS_CMD_MICRO = """
    preprocess \
    -infile #VEP# \
    -ds #DS_FILE# \
    -out #OUT# \
    -vep2mti
"""

PREPROCESS_CMD_FULL = """
    preprocess \
    -infile #VEP# \
    -out #OUT# \
    -vep2mti
"""

MICROANNOT_CMD = """
    annot \
    -vcf #VCF# \
    -ds  #MICRO_DS_FILE# \
    -out #OUT# \
    -sourcefile #MICRO_SOURCE_FILE# #MICROANNOT_VEP_MTI#\
    -split_multi_allelic_variant \
    -genoinfo \
    -use_raw_source
"""

FULLANNOT_CMD = """
    annot \
    -vcf #VCF# \
    -ds  #FULL_DS_FILE# \
    -out #OUT# \
    -sourcefile #FULL_SOURCE_FILE# #FULLANNOT_VEP_MTI#\
    -hg19 \
    -chain #CHAINFILE# \
    -clean_tag MUTANNO SpliceAI CLINVAR gnomADgenome
"""

TEST_VERSION = conf.getNow()

MICROANNOT_VEP_MTI = "vep_output_microannot.mti.gz"
FULLANNOT_VEP_MTI = "vep_output_fullannot.mti.gz"

################################
## MICRO ANNOTATION
################################

def test_preprocess_vep_to_mti_micro():
    cmd = PREPROCESS_CMD_MICRO
    opt = {}
    opt['VEP'] = "#TESTDATA_PATH#/vep_output.vcf.gz"
    opt['DS_FILE'] = "#MICRO_DS_FILE#"
    opt['OUT'] = "#TESTOUT_PATH#/" + MICROANNOT_VEP_MTI.replace('.mti.gz', '.mti')
    opt['PREV_OUT'] = "#TESTDATA_PATH#/" + MICROANNOT_VEP_MTI
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
    mutanno.cli()

    conf.tabixgz(opt['OUT'])
    conf.comp_previous_out(opt['OUT'], opt['PREV_OUT'])


def check_vep(annot_vcf, source_mti):
    m = {}
    c = {}
    for line in mutanno.util.file_util.gzopen(source_mti):
        line = mutanno.util.file_util.decodeb(line)
        if line[0] != '#':
            arr = line.split('\t')
            if 'VEP=' in arr[7]:
                varkey = arr[0] + '_' + arr[1] + '_' + arr[3] + '_' + arr[4]
                m[varkey] = arr[7].strip()
                c[varkey] = 1

    no_match = 0
    for line in mutanno.util.file_util.gzopen(annot_vcf):
        line = mutanno.util.file_util.decodeb(line)
        if line[0] != '#':
            arr = line.split('\t')
            varkey = arr[0].replace('chr','') + '_' + arr[1] + '_' + arr[3] + '_' + arr[4]
            try:
                m[varkey]
                # print('>',varkey)
                for s1 in arr[7].split(';'):
                    if 'VEP=' in s1:
                        no_match += 1
                        c[varkey] += 1
            except KeyError:
                pass
    # print(no_match, len(m.keys()))
    # for v1 in c.keys():
    #     print(v1, c[v1])
    assert no_match == len(m.keys())
    

def test_microannot_with_multiple_mti():
    cmd = MICROANNOT_CMD
    opt = {}
    opt['VCF'] = '#TESTDATA_PATH#/TRIO_GAPFIS8ZSPEO.target.add_novep.vcf.gz'
    opt['PREV_OUT'] = '#TESTDATA_PATH#/TRIO_GAPFIS8ZSPEO.target.add_novep.vcf.gz.microannot_v0.4.4.vcf.gz'
    opt['OUT'] = '#TESTOUT_PATH#/TRIO_GAPFIS8ZSPEO.target.add_novep.vcf.gz.microannot_v0.4.4.'+TEST_VERSION+'.vcf'
    opt['MICROANNOT_VEP_MTI'] = "#TESTOUT_PATH#/" + MICROANNOT_VEP_MTI
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
    check_vep(opt['OUT'], opt['MICROANNOT_VEP_MTI'])
    # conf.check_annotation_field(opt['OUT'], ['VEP', 'SAMPLEGENO'])
    conf.comp_previous_out(opt['OUT'], opt['PREV_OUT'])
    # conf.validate_annotvcf(opt['OUT'], opt['MICRO_DS_FILE'], sys.argv)


################################
## FULL ANNOTATION 
################################

def test_preprocess_vep_to_mti_full():
    cmd = PREPROCESS_CMD_FULL
    opt = {}
    opt['VEP'] = "#TESTDATA_PATH#/vep_output.vcf.gz"
    opt['OUT'] = "#TESTOUT_PATH#/" + FULLANNOT_VEP_MTI.replace('.mti.gz', '.mti')
    opt['PREV_OUT'] = "#TESTDATA_PATH#/" + FULLANNOT_VEP_MTI
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
    mutanno.cli()

    conf.tabixgz(opt['OUT'])
    conf.comp_previous_out(opt['OUT'], opt['PREV_OUT'])

def test_fullannot_with_multiple_mti():
    cmd = FULLANNOT_CMD
    opt = {}
    # opt['VCF'] = '#TESTOUT_PATH#/TRIO_GAPFIS8ZSPEO.target.add_novep.vcf.gz.microannot_v0.4.4.'+TEST_VERSION+'.vcf'
    opt['VCF'] = "#TESTDATA_PATH#/TRIO_GAPFIS8ZSPEO.target.add_novep.vcf.gz.microannot_v0.4.4.vcf.gz"
    # opt['VCF'] = "#TESTDATA_PATH#/TRIO_GAPFIS8ZSPEO.target.add_novep.vcf.gz.microannot_v0.4.4-1.vcf"
    # opt['VCF'] = "#TESTDATA_PATH#/TRIO_GAPFIS8ZSPEO.target.add_novep.vcf.gz.microannot_v0.4.4-2.vcf"
    
    opt['PREV_OUT'] = '#TESTDATA_PATH#/TRIO_GAPFIS8ZSPEO.target.add_novep.vcf.gz.microannot_v0.4.4.fullannot_v0.4.8.vcf'
    opt['OUT'] = '#TESTOUT_PATH#/TRIO_GAPFIS8ZSPEO.target.add_novep.vcf.gz.microannot_v0.4.4.fullannot_v0.4.8.'+TEST_VERSION+'.vcf'
    opt['FULLANNOT_VEP_MTI'] = "#TESTOUT_PATH#/" + FULLANNOT_VEP_MTI
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
    check_vep(opt['OUT'], opt['FULLANNOT_VEP_MTI'])

    stat = conf.validate_annotvcf(opt['OUT'], opt['FULL_DS_FILE'], sys.argv)
    prev_stat = conf.validate_annotvcf(opt['PREV_OUT'], opt['FULL_DS_FILE'], sys.argv)
    # print("STAT:", stat)
    # print("PREV STAT:", prev_stat)
    print ('Key', 'Stat', 'Prev stat')
    for k1 in prev_stat.keys():
        print(k1, stat[k1], prev_stat[k1])
        assert stat[k1] == prev_stat[k1]

    conf.comp_previous_out(opt['OUT'], opt['PREV_OUT'])




if __name__ == "__main__":
    test_preprocess_vep_to_mti_micro()
    test_microannot_with_multiple_mti()

    test_preprocess_vep_to_mti_full()
    test_fullannot_with_multiple_mti()


#
