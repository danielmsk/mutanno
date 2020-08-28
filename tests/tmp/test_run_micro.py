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


def comp_previous_out(out, prevout):
    outcont = file_util.fileOpen(out)
    prevoutcont = file_util.fileOpen(prevout)
    assert outcont == prevoutcont, "Not matched: " + out + " " + prevout

def check_vcf_validator(out):
    # assert test_conf.run_vcf_validator(out) == test_conf.DEFAULT_VCF_VALIDATOR_MSG
    assert test_conf.run_vcf_validator(out) == ''

def validate_annotvcf(vcf, dsjson, opt):
    va = AnnotVCFValidator(vcf)
    va.set_datastructure(jsonfile=dsjson)
    va.set_option(opt)
    va.validate()

def validate_annottsi(tsi, dsjson, opt):
    va = AnnotTSIValidator(tsi)
    va.set_datastructure(jsonfile=dsjson)
    va.set_option(opt)
    va.validate()


def load_googlesheet(googlesheet):
    rst = {}
    i = 0
    for line in open(googlesheet):
        if len(line.strip()) > 0:
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            if line.strip()[0] != "#":
                if i == 0 :
                    header = arr
                else:
                    d = {}
                    for i in range(len(header)):
                        header[i] = header[i].strip()
                        if header[i] != '':
                            d[header[i]] = arr[i].strip()
                    rst[d['field_name']] = d
                i += 1
    return rst

def load_dsjson(dsjson):
    js = json.load(open(dsjson,'r'))
    rst = {}
    for s1 in js['source']:
        # print(s1['name'])
        if dv(s1, 'is_available', True):
            for f1 in s1['fields']:
                fname = dv(f1, 'name2', f1['name'])
                k1 = s1['name'].lower() + '_' + fname.lower()
                f1['subembedded'] = dv(s1, 'subembedded', '')
                rst[k1] = f1
    return rst


def test_microannot_run():
    ds = test_conf.get_ds()
    # for n1 in [10, 100, 1000]:
    # for n1 in [10, 100]:
    # for n1 in [100]:
    for n1 in [10]:
        for r1 in range(1,6):
            vcf = "data/test_trio_"+str(n1)+"_"+str(r1)+".vcf.gz"
            # vcf = "data/test_novocaller_"+str(n1)+"_"+str(r1)+".vcf"
            out = "out/" + vcf.split('/')[-1] + '.microannot.vcf'
            prevout = "data/" + out.split('/')[-1]
            # dsjson = ds['microannot']['ds']
            dsjson = ds['microannot']['clean']

            arg = ['mutanno']
            arg.append('annot')
            arg.extend(['-vcf',vcf])
            arg.extend(['-ds',dsjson])
            arg.extend(['-out',out])
            arg.extend(['-sourcefile',test_conf.SOURCEFILE['microannot']])
            arg.append('-split_multi_allelic_variant')
            arg.append('-single_source_mode')
            # arg.extend(['-genoinfo', 'NA12877_sample', 'NA12878_sample', 'NA12879_sample'])
            arg.extend(['-genoinfo'])

            sys.argv = arg
            print('>>command:',' '.join(sys.argv))
            mutanno.cli()

            check_vcf_validator(out)
            comp_previous_out(out, prevout)
            validate_annotvcf(out, dsjson, arg)
            break
            

def test_novocaller_vcf_format():
    # for n1 in [10, 100, 1000]:
    for n1 in [10]:
        for r1 in range(1,6):
            out = "data/test_novocaller_"+str(n1)+"_"+str(r1)+".vcf"
            # check_vcf_validator(out)
            # validate_annotvcf(out)
            break


def print_command():
    print('>>command: ' + ' '.join(sys.argv))


if __name__ == "__main__":
    
    test_microannot_run()
    



# 