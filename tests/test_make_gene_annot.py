import sys
import shlex
import filecmp
import json
import time
import conf

import mutanno

prog = "mutanno"

MK_GENEANNOT_CMD = """
    makedata \
    -ds  #GENE_DS_FILE# \
    -out #OUT# \
    -vartype CODING_GENE_MAIN_CHROM \
    -outtype json
"""

TEST_VERSION = conf.getNow()

def test_gene_annot():
    cmd = MK_GENEANNOT_CMD
    opt = {}
    opt['OUT'] = '#TESTOUT_PATH#/mvp_gene_datasource_v0.4.6.coding_gene_main_chrom.'+TEST_VERSION
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

    conf.tabixgz(opt['OUT'])
    conf.check_vcf_validator(opt['OUT']+'.gz')

if __name__ == "__main__":
    test_gene_annot()
    
