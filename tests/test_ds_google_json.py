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
    skiplist = ["VCF","SAMPLEGENO","MULTIALLELE","MUTANNO","HG19","VEP", "COMHET","HGVS"]
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
    


def print_command():
    print('>>command: ' + ' '.join(sys.argv))


if __name__ == "__main__":
    test_sync_googlesheet_dsjson()
    
    
    



