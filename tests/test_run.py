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

def check_fields_property_googlesheet_dsjson(g1, d1, k1, is_exception=True):
    if (g1[k1]['do_import'] == "Y") != dv (d1[k1], 'is_available', True):
        msg = k1 + " field's do_import/is_available is not matched."
        if is_exception:
            raise Exception(msg)
        else:
            print("WARN:"+msg)

    if g1[k1]['field_type'] != d1[k1]['type']:
        msg = k1 + " field's type is not matched."
        if is_exception:
            raise Exception(msg)
        else:
            print("WARN:"+msg)           

    if (g1[k1]['is_list'] == "Y") != dv(d1[k1], 'is_list', False):
        msg = k1 + " field's is_list is not matched."
        if is_exception:
            raise Exception(msg)
        else:
            print("WARN:"+msg)

    if (g1[k1]['sub_embedding_group'].strip() != ""):
        e1 = json.loads(g1[k1]['sub_embedding_group'].strip())
        if e1['key'] != dv(d1[k1], 'subembedded', ''):
            msg = k1 + " field's subembedded is not matched."
            if is_exception:
                raise Exception(msg)
            else:
                print("WARN:"+msg)

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


def test_sync_googlesheet_dsjson():
    ds = test_conf.get_ds()
    dsjson = ds['fullannot']['ds']
    dsjson_micro = ds['microannot']['ds']
    googlesheet = dsjson.replace(".json","") + ".googlesheet.txt"
    g1 = load_googlesheet(googlesheet)
    d1 = load_dsjson(dsjson)
    m1 = load_dsjson(dsjson_micro)

    i = 0
    print("DSJSON(full) -> Google")
    skiplist = ["ENSEMBLANNOT"]
    for k1 in d1.keys():
        if dv (d1[k1], 'is_available', True) and k1.split('_')[0].upper() not in skiplist:
            i += 1
            # print('\t',i, k1, '..')
            try:
                g1[k1]
            except KeyError:
                raise Exception(k1 + " field is not in googlesheet.")
            check_fields_property_googlesheet_dsjson(g1, d1, k1)

    print("\t", i, " dsjson(fullannot) fields checked...")

    i = 0
    print("Google -> DSJSON(full)")
    skiplist = ["VCF","SAMPLEGENO","MULTIALLELE","MUTANNO","HG19","VEP", "COMHET"]
    for k1 in g1.keys():
        if (g1[k1]['do_import'] == "Y") and (g1[k1]['source_name'] not in skiplist) and ( g1[k1]['embedded_field'] != "Y"):
            i += 1
            # print('\t',i, k1)
            try:
                d1[k1]
            except KeyError:
                raise Exception(k1 + " field is not in DS json.")

            check_fields_property_googlesheet_dsjson(g1, d1, k1)
    print("\t", i, " google(fullannot) fields checked...")

    i = 0
    print("DSJSON(micro) -> Google")
    for k1 in m1.keys():
        if dv (m1[k1], 'is_available', True) and k1.startswith('vep_') :
            i += 1
            # print('\t',i, k1, "..")
            try:
                g1[k1]
            except KeyError:
                raise Exception(k1 + " field is not in googlesheet.")


            check_fields_property_googlesheet_dsjson(g1, m1, k1, False)
    print("\t", i, " dsjson(microannot) fields checked...")

    # i = 0
    # print("Google -> DSJSON(micro)")
    # for k1 in g1.keys():
    #     if (g1[k1]['do_import'] == "Y") and (g1[k1]['source_name'] == "VEP") and ( g1[k1]['embedded_field'] != "Y"):
    #         i += 1
    #         try:
    #             m1[k1]
    #         except KeyError:
    #             raise Exception(k1 + " field is not in DS json.")
    #         check_fields_property_googlesheet_dsjson(g1, m1, k1)
    # print("\t", i, " google(microannot, VEP) fields checked...")
    

def test_microannot_run():
    ds = test_conf.get_ds()
    for n1 in [10, 100, 1000]:
    # for n1 in [10, 100]:
    # for n1 in [10]:
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
            # arg.append('-single_source_mode')
            # arg.extend(['-genoinfo', 'NA12877_sample', 'NA12878_sample', 'NA12879_sample'])
            arg.extend(['-genoinfo'])

            sys.argv = arg
            print('>>command:',' '.join(sys.argv))
            mutanno.cli()

            check_vcf_validator(out)
            # comp_previous_out(out, prevout)
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


def test_fullannot_run():
    ds = test_conf.get_ds()
    for n1 in [10]:
        for r1 in range(1,6):
            # vcf = "data/test_novocaller_"+str(n1)+"_"+str(r1)+".vcf"
            # vcf = "data/test_novocaller_"+str(n1)+"_"+str(r1)+".vcf" + '.microannot.vcf'
            vcf = "out/test_trio_"+str(n1)+"_"+str(r1)+".vcf.gz" + '.microannot.vcf'
            print(vcf)
            out = "out/" + vcf.split('/')[-1] + '.fullannot.vcf'
            prevout = "data/" + out.split('/')[-1]
            dsjson = ds['fullannot']['ds']

            arg = ['mutanno']
            arg.append('annot')
            arg.extend(['-vcf',vcf])
            arg.extend(['-ds',dsjson])
            arg.extend(['-out',out])
            arg.extend(['-outtype', 'vcf', 'json'])
            arg.append('-hgvs')
            arg.append('-variant_class')
            arg.append('-hg19')
            arg.extend(['-chain',test_conf.CHAINFILE])
            arg.extend(['-clean_tag','SpliceAI','CLINVAR','gnomADgenome'])
            sys.argv = arg

            # sys.argv.extend(['-sourcefile',test_conf.SOURCEFILE['fullannot']])
            # print('>>command:',' '.join(sys.argv))
            # mutanno.cli()
            # check_vcf_validator(out)
            # comp_previous_out(out, prevout)
            # validate_annotvcf(out, dsjson,arg)
            break


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


def print_command():
    print('>>command: ' + ' '.join(sys.argv))

def test_make_microannot_sourcefile():
    ds = test_conf.get_ds()
    # dsjson = ds['fullannot']['ds']
    key = "microannot"
    dsjson = ds[key]['ds']
    out = "out/"+key+".datasource_test1.tsi"

    arg = ['mutanno']
    arg.append('makedata')
    arg.extend(['-ds',dsjson])
    arg.extend(['-out',out])
    arg.extend(['-vartype','SNV'])
    # arg.extend(['-region','1:22500-22501']) # gnomAD only variant
    # arg.extend(['-region','1:952029-952609']) # clinvar variant
    # arg.extend(['-region','1:952029-952030']) # clinvar variant
    # arg.extend(['-region','1:10163-10174']) # clinvar variant
    # arg.extend(['-region','1:15034-15035']) # clinvar varian
    arg.extend(['-region','1:10001-11000']) # clinvar variant
    # arg.extend(['-region','5:109000001-109000002']) 
    
    # sys.argv = arg
    # print_command()
    # mutanno.cli()

    # validate_annottsi(out, dsjson, arg)

def test_make_fullannot_sourcefile():
    ds = test_conf.get_ds()
    # dsjson = ds['fullannot']['ds']
    key = "fullannot"
    dsjson = ds[key]['ds']
    out = "out/"+key+".datasource_test1.tsi"

    arg = ['mutanno']
    arg.append('makedata')
    arg.extend(['-ds',dsjson])
    arg.extend(['-out',out])
    arg.extend(['-vartype','SNV'])
    arg.extend(['-blocksize','5000'])
    # arg.extend(['-region','1:10001-10002']) # clinvar variant
    # arg.extend(['-region','1:22500-22501']) # gnomAD only variant
    # arg.extend(['-region','1:952029-952659']) # clinvar variant
    # arg.extend(['-region','1:952029-952030']) # clinvar variant
    arg.extend(['-region','1:10163-20174']) # clinvar variant
    # arg.extend(['-region','1:15034-15035']) # clinvar variant
    # arg.extend(['-region','1:39984-40000']) # clinvar variant
    
    # arg.extend(['-region','1:59900001-59900010']) 
    
    # sys.argv = arg
    # print_command()
    # mutanno.cli()

    # validate_annottsi(out, dsjson, arg)


if __name__ == "__main__":
    test_microannot_run()
    # test_novocaller_vcf_format()
    # test_sync_googlesheet_dsjson()
    # test_fullannot_run()
    # test_allannot_run()
    # test_make_microannot_sourcefile()
    # test_make_fullannot_sourcefile()
    



