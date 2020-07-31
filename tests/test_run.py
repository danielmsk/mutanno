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
            # print('>>command:',' '.join(sys.argv))
            # mutanno.cli()
            # check_vcf_validator(out)
            # comp_previous_out(out, prevout)
            # validate_annotvcf(out, dsjson, arg)
            
            

def test_novocaller_vcf_format():
    # for n1 in [10, 100, 1000]:
    for n1 in [10]:
        for r1 in range(1,6):
            out = "data/test_novocaller_"+str(n1)+"_"+str(r1)+".vcf"
            # check_vcf_validator(out)
            # validate_annotvcf(out)
            break


'''
mutanno annot \
    -vcf test.vcf
    -out test.annot.vcf
    -hgvs
    -variant_class
    -hg19
    -clean_tag VEP SpliceAI CLINVAR gnomADgenome
    -signle_source_mode
    -chain hg38tohg19.chain
    -ds datastructure_fullannot_v0.4.6.json
    -sourcefile fullannot.source.tsi.gz
'''

def test_fullannot_run():
    ds = test_conf.get_ds()
    for n1 in [10]:
        for r1 in range(1,6):
            # vcf = "data/test_novocaller_"+str(n1)+"_"+str(r1)+".vcf"
            # vcf = "data/test_novocaller_"+str(n1)+"_"+str(r1)+".vcf" + '.microannot.vcf'
            # vcf = "out/test_trio_"+str(n1)+"_"+str(r1)+".vcf.gz" + '.microannot.vcf'
            # vcf = "/home/mk446/mutanno/TEST/0616/GAPFIP83PL7E_test1.vcf.gz"
            # vcf = "/home/mk446/mutanno/TEST/0616/GAPFIP83PL7E_test2.vcf.gz"
            vcf = "/home/mk446/mutanno/TEST/0616/GAPFI3JX5D2J.vcf.gz"
            # vcf = "/home/mk446/mutanno/TEST/0616/GAPFI3JX5D2J_test1.vcf"
            
            # vcf = "/home/mk446/mutanno/TEST/0616/GAPFI3JX5D2J_chr9.vcf.gz"
            # vcf = "/home/mk446/mutanno/TEST/0616/GAPFI3JX5D2J_test1.vcf.gz"
            # vcf = "/home/mk446/mutanno/TEST/0616/GAPFI3JX5D2J_test2.vcf.gz"
            # vcf = "/home/mk446/mutanno/TEST/0616/GAPFI3JX5D2J_test3.vcf.gz"
            # vcf = "/home/mk446/mutanno/TEST/0616/GAPFI3JX5D2J_test4.vcf.gz"
            # out = vcf + '.fullannot.vcf'
            out = vcf + '.fullannot.wbfilter.vcf'
            # out = "/home/mk446/mutanno/TEST/0616/GAPFIRHN9YOZ.vcf"
            # out = "out/" + vcf.split('/')[-1] + '.fullannot.vcf'
            prevout = "data/" + out.split('/')[-1]
            dsjson = ds['fullannot']['clean']
            # dsjson = "/home/mk446/mutanno/SRC/tests/data/archive/datastructure_fullannot_v0.4.6.json"
            print("DS:",dsjson)
            # sourcefile = "/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/v0.4.6_wbfilter_split/BGZIP/fullannot_datasource.09_06.v0.4.6_200617.tsi.gz"
            # sourcefile = "/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/v0.4.6_wbfilter_split/BGZIP/22_00.tsi.gz"
            # sourcefile = "/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/v0.4.6_wbfilter_split/BGZIPLINK/#CHROM#_00.tsi.gz"
            sourcefile = "/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/v0.4.6_wbfilter_split/merged.mti.gz"
            # sourcefile = "/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/GAPFI3JX5D2J_fullannot_source.sorted.rmdup.tsi.gz"
            # sourcefile = '/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/v0.4.6_wbfilter/fullannot_datasource.#CHROM#.v0.4.6_200617.tsi.gz'
            # sourcefile = '/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/v0.4.6_all/BGZIP/fullannot_datasource.#CHROM#.v0.4.6_200617.tsi.gz'
            # sourcefile = '/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/v0.4.6_wbfilter/GAPFIR6C5TC7.mti.gz'
            # sourcefile = '/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/v0.4.6_all/x.chrom.tsi.gz'

            arg = ['mutanno']
            arg.append('annot')
            arg.extend(['-vcf',vcf])
            arg.extend(['-ds',dsjson])
            arg.extend(['-out',out])
            arg.extend(['-outtype', 'vcf', 'json'])
            # arg.extend(['-outtype', 'json'])
            # arg.extend(['-outtype', 'vcf'])
            arg.extend(['-sourcefile',sourcefile])
            arg.append('-hg19')
            arg.extend(['-chain',test_conf.CHAINFILE])

            # arg.append('-hgvs')
            # arg.append('-variant_class')
            # arg.extend(['-clean_tag','VEP','SpliceAI','CLINVAR','gnomADgenome'])

            arg.extend(['-clean_tag','MUTANNO','SpliceAI','CLINVAR','gnomADgenome'])

            # arg.append('-single_source_mode')
            # arg.extend(['-sourcefile', test_conf.SOURCEFILE['fullannot']])

            sys.argv = arg
            
            print('>>command:',' '.join(sys.argv))
            mutanno.cli()
            check_vcf_validator(out)
            validate_annotvcf(out, dsjson,arg)
            # comp_previous_out(out, prevout)
            break
        break


'''
def test_allannot_run():
    ds = test_conf.get_ds()
    for n1 in [10]:
        for r1 in range(1,6):
            vcf = "data/test_trio_"+str(n1)+"_"+str(r1)+".vcf"
            out = "out/" + vcf.split('/')[-1] + '.allannot.vcf'
            prevout = "data/" + out.split('/')[-1]
            dsjson = ds['allannot']['ds']

            arg = ['mutanno']
            arg.append('annot')
            arg.extend(['-vcf',vcf])
            arg.extend(['-ds',dsjson])
            arg.extend(['-out',out])
            arg.extend(['-outtype', 'vcf', 'json'])
            arg.append('-hgvs')
            arg.append('-hg19')
            arg.append('-variant_class')
            arg.append('-split_multi_allelic_variant')
            arg.extend(['-genoinfo', 'NA12877_sample', 'NA12878_sample', 'NA12879_sample'])
            arg.extend(['-chain',test_conf.CHAINFILE])
            # arg.extend(['-clean_tag','SpliceAI','CLINVAR','gnomADgenome'])
            # sys.argv.extend(['-sourcefile',test_conf.SOURCEFILE['fullannot']])
            # print('>>command:',' '.join(sys.argv))
            sys.argv = arg
            # mutanno.cli()
            # check_vcf_validator(out)
            # comp_previous_out(out, prevout)
            # validate_annotvcf(out, dsjson, arg)
            # break
'''

def print_command():
    print('>>command: ' + ' '.join(sys.argv))


if __name__ == "__main__":
    
    # test_novocaller_vcf_format()
    # test_microannot_run()
    test_fullannot_run()
    # test_allannot_run()
    



# 