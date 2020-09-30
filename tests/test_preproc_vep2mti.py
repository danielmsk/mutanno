import sys
import conf
import shlex

sys.path.append('..')
from src import mutanno

prog = "mutanno"


PREPROCESS_CMD_FULL = """
    preprocess \
    -infile #VEP# \
    -out #OUT# \
    -vep2mti
"""

FULLANNOT_VEP_MTI = "vep_output2.fullannot.mti.gz"


def check_vep_field(vep2mti_file):
    for line in mutanno.util.file_util.gzopen(vep2mti_file):
        line = mutanno.util.file_util.decodeb(line)
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        varkey = arr[0] + ':' + arr[1] + '_' + arr[3] + '>' + arr[4]
        if line[0] == '#':
            fieldnames = arr[-1].replace('VEP=', '').split('|')
        else:
            if 'VEP=' in arr[-1]:
                for section in arr[-1].split(','):
                    arr2 = section.split('|')
                    assert len(fieldnames) == len(arr2), "The number of fields is not matched. " + \
                        str(len(fieldnames)) + ' != ' + str(len(arr2)) + " in " + varkey
            else:
                raise Exception("no VEP field in " + varkey)


def test_preproc_vep2mti_minus_deletion():
    cmd = PREPROCESS_CMD_FULL
    opt = {}
    # opt['VEP'] = "#TESTDATA_PATH#/vep_output2.vcf.gz"
    # opt['VEP'] = "#TESTDATA_PATH#/multiallelic_star_deletion.vep.vcf.gz"
    opt['VEP'] = "#TESTDATA_PATH#/multiallelic_minus_deletion.vep.vcf.gz"
    opt['OUT'] = "#TESTOUT_PATH#/" + opt['VEP'].split('/')[-1].replace('.vcf.gz', '.mti')
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
    check_vep_field(opt['OUT'])
    
    # conf.comp_previous_out(opt['OUT'], opt['PREV_OUT'])


def test_preproc_vep2mti_star_deletion_missing_vep_annotation_for_star():
    cmd = PREPROCESS_CMD_FULL
    opt = {}
    # opt['VEP'] = "#TESTDATA_PATH#/vep_output2.vcf.gz"
    opt['VEP'] = "#TESTDATA_PATH#/multiallelic_star_deletion.vep.vcf.gz"
    # opt['VEP'] = "#TESTDATA_PATH#/multiallelic_minus_deletion.vep.vcf.gz"
    opt['OUT'] = "#TESTOUT_PATH#/" + opt['VEP'].split('/')[-1].replace('.vcf.gz', '.mti')
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
    check_vep_field(opt['OUT'])

if __name__ == "__main__":
    test_preproc_vep2mti_minus_deletion()
    test_preproc_vep2mti_star_deletion_missing_vep_annotation_for_star()

