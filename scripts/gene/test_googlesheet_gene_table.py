#!/usr/bin/env python
# -*- coding: utf-8 -*-
# check_googlesheet_gene_table.py
# made by Daniel Minseok Kwon
# 2020-04-27 10:56:56
#########################
import sys
import os
import json
import re
import pytest
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path = "/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path = "/ms1/bin/python_lib"
else:
    sys_path = "/home/mk446/bin/python_lib"
sys.path.append(sys_path)

CURRENTPATH = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(CURRENTPATH, "../")))

import proc_util
import file_util
import preproc_util

GS_FIELD_NAME = {}
GS_FIELD_NAME['FIELD_TYPE'] = 'field_type (string, integer, boolean, number)'
GS_FIELD_NAME['IS_LIST'] = 'is_list'
GS_FIELD_NAME['DO_IMPORT'] = 'do_import'
GS_FIELD_NAME['ENUM_LIST'] = 'enum_list (limit values to this list)'
GS_FIELD_NAME['PATTERN'] = 'pattern'

datapath = preproc_util.DATASOURCEPATH + '/MUTANOANNOT/'
VERSION = "0.4.4"


def load_googlesheet_ds(googlesheet_ds):
    h = {}
    header = []
    gs = {}
    subembed = {}
    for line in open(googlesheet_ds):
        if line.strip()[0] != '#':
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            if arr[0] == "no":
                for k in range(len(arr)):
                    h[arr[k]] = k
                    header.append(arr[k].strip())
            else:
                field_name = arr[h['Field_Name']].strip()

                assert (field_name in gs.keys()) == False, "Duplicate field: " + field_name
                # try:
                #     gs[field_name]
                #     print('Duplicate Error: ' + field_name)
                # except KeyError:
                gs[field_name] = {}


                d = {}
                for k in range(len(arr)):
                    d[header[k]] = arr[k].strip()

                subembed_key = arr[h["sub_embedding_group (SUBEMBEDDED OBJECT GROUPING)"]].strip()
                if subembed_key != "":
                    try:
                        subembed[subembed_key]
                    except KeyError:
                        subembed[subembed_key] = {}
                    subembed[subembed_key][field_name] = d
                else:
                    gs[field_name] = d
    return gs, subembed



class GeneAnnotData:
    def __init__(self, gene_annot_data_file):
        self.fp = file_util.gzopen(gene_annot_data_file)
        self.init_flag = False
        self.depth = 0

    def __iter__(self):
        return self

    def __next__(self):
        cont = ""
        i = 0
        while True:
            i += 1
            if i > 1000000:
                break
            w = file_util.decodeb(self.fp.read(1))

            if w == "{":
                self.depth += 1

            if w == "}":
                self.depth -= 1
                if self.depth == 0:
                    cont += w
                    break
            if self.depth > 0:
                cont += w
        if cont.strip() == "":
            raise StopIteration
        return json.loads(cont)

def check_field_type(logkey, gsdata, rec_value):
    global GS_FIELD_NAME
    
    assert gsdata[GS_FIELD_NAME['IS_LIST']] in ['Y','N'], logkey

    if gsdata[GS_FIELD_NAME['IS_LIST']] == 'Y':
        rec_list = rec_value
    elif gsdata[GS_FIELD_NAME['IS_LIST']] == 'N':
        rec_list = [rec_value]

    for value1 in rec_list:
        log = logkey
        field_type = gsdata[GS_FIELD_NAME['FIELD_TYPE']]
        gs_field_type = ''
        if field_type == "string":
            gs_field_type = str
        if field_type == "integer":
            gs_field_type = int
        if field_type == "number":
            gs_field_type = float
        if field_type == "boolean":
            gs_field_type = bool
        
        assert type(value1) == gs_field_type, logkey

def check_do_import(logkey, gsdata):
    global GS_FIELD_NAME
    assert gsdata[GS_FIELD_NAME['DO_IMPORT']] == 'Y', logkey 

def check_enum_list(logkey, gsdata, rec_value):
    global GS_FIELD_NAME
    enum_list_cont = gsdata[GS_FIELD_NAME['ENUM_LIST']].strip()
    if enum_list_cont != "":
        enum_list = []
        for el in enum_list_cont.split(','):
            enum_list.append(el.strip())
        assert rec_value in enum_list, logkey

def check_pattern(logkey, gsdata, rec_value):
    global GS_FIELD_NAME

    if gsdata[GS_FIELD_NAME['IS_LIST']] == 'Y':
        rec_list = rec_value
    elif gsdata[GS_FIELD_NAME['IS_LIST']] == 'N':
        rec_list = [rec_value]

    for value1 in rec_list:
        pattern = gsdata[GS_FIELD_NAME['PATTERN']].strip()
        if pattern != "":
            pattern = pattern.replace('\\\\','\\')
            rst = re.match(pattern, value1)
            if pattern[-1] == "$":
                assert rst.group(0) == value1, logkey
            else:
                assert rst.group(0) in value1, logkey

@pytest.mark.parametrize("version",[
    (VERSION)
    ])

def test_googlesheet_gene_table(version):
    googlesheet_ds = preproc_util.TESTPATH + "/data/datastructure_gene_v"+version+"ds.googlesheet.txt"
    genetable_data = datapath + "mvp_gene_datasource_v"+version+".coding_gene_main_chrom.json.gz"

    gs, subembed = load_googlesheet_ds(googlesheet_ds)

    gadata = GeneAnnotData(genetable_data)
    i = 0
    for rec in gadata:
        i += 1

        for field_name in rec.keys():
            # print('['+str(i)+']' + rec['gene_symbol']+':'+field_name)
            if field_name in gs.keys():
                logkey = '['+str(i)+']' + rec['gene_symbol']+':'+field_name

                check_field_type(logkey, gs[field_name], rec[field_name])
                check_do_import(logkey, gs[field_name])
                check_enum_list(logkey, gs[field_name], rec[field_name])
                check_pattern(logkey, gs[field_name], rec[field_name])
                
            else:
                for subfield in rec[field_name]:
                    for subfieldname in subfield.keys():
                        logkey = '['+str(i)+']' + rec['gene_symbol']+':'+field_name+':'+subfieldname
                        check_field_type(logkey, subembed[field_name][subfieldname], subfield[subfieldname])
                        check_do_import(logkey, subembed[field_name][subfieldname])
                        check_enum_list(logkey, subembed[field_name][subfieldname], subfield[subfieldname])
                        check_pattern(logkey, subembed[field_name][subfieldname], subfield[subfieldname])

        if i > 200:
            pass


if __name__ == "__main__":
    test_googlesheet_gene_table(version)

    
    



