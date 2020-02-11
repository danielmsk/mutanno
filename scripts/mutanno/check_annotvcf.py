#!/usr/bin/env python
# -*- coding: utf-8 -*-
# check_annotvcf.py
# made by Daniel Minseok Kwon
# 2020-02-06 16:22:07
#########################
import sys
import os
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path = "/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path = "/ms1/bin/python_lib"
else:
    sys_path = "/home/mk446/bin/python_lib"
sys.path.append(sys_path)


def get_source_list(mapping_table_col, mapping_table):
    source_list = []
    source_version_map = {}
    source_fields = {}
    for i in range(len(mapping_table_col['SOURCE_NAME'])):
        fieldname = mapping_table_col['VCF NAME (field name on ann vcf)'][i]

        if mapping_table_col['SOURCE_NAME'][i] not in source_list:
            source_name = mapping_table_col['SOURCE_NAME'][i]
            source_list.append(source_name)
            source_version_map[source_name] = mapping_table_col['SOURCE VERSION'][i]

        try:
            tmp = source_fields[source_name]
        except KeyError:
            source_fields[source_name] = {}

        if fieldname not in source_fields[source_name].keys():
            if mapping_table_col['SCOPE (sample_variant/variant/gene)'][i] == "variant":
                source_fields[source_name][fieldname] = mapping_table[fieldname]
        
    return (source_list, source_version_map, source_fields)


def check_source_mutanno_header(line, source_list, source_version_map):
    arr = line.split(',')
    sname = arr[0].replace("##MUTANNO=<ID=","")
    sversion = arr[1].replace("Version=","").replace('"','')
    # print(line)
    if sname not in source_list:
        if sname not in IGNORESOURCE:
            print ("ERROR: only in VCF:", sname)
    else:
        if source_version_map[sname] != sversion:
            print ("ERROR: different source version:", sname, sversion + " (VCF) != " + source_version_map[sname] + " (mapping table)")
    return sname

def check_source_fields_header(line, source_list, source_fields):
    if ' Format:' in line:
        sname = line.split(',')[0].replace("##INFO=<ID=","").strip()
        vcf_fields = line.split(' Format:')[-1].replace("'","").replace('">','').strip().split('|')

        if sname in source_list:
            doc_fields = list(source_fields[sname].keys())
            if len(vcf_fields) != len(doc_fields):
                print ("ERROR: not matching field number of " + sname + ": " + str(len(vcf_fields)) + " vs. " + str(len(doc_fields)))

            vcf_fields2 = []
            for vcf_field in vcf_fields:
                vcf_field2 = sname.lower() + '_' + vcf_field.replace(' ','_').lower()
                vcf_fields2.append(vcf_field2)
                if vcf_field2 not in doc_fields:
                    print("ERROR: field is not in mapping table: " + vcf_field2)

            for doc_field in doc_fields:
                if doc_field not in vcf_fields2:
                    print("ERROR: field is not in VCF: " + doc_field)


def check_annotvcf_header(annotvcf, source_list, source_version_map, source_fields):
    vcf_source_list = []
    for line in file_util.gzopen(annotvcf):
        if annotvcf.endswith('.gz'):
            line = line.decode('UTF-8')
        if line[0] == "#":
            if line[:len('##MUTANNO=<ID=')] == "##MUTANNO=<ID=":
                vcf_source_list.append(check_source_mutanno_header(line, source_list, source_version_map))
            if line[:len('##INFO=<ID=')] == "##INFO=<ID=":
                # print(line.strip())
                check_source_fields_header(line, source_list, source_fields)
        else:
            break

    for sname in source_list:
        if sname not in vcf_source_list:
            if sname not in IGNORESOURCE:
                print ("ERROR: only in mapping table:", sname)


def check_annotvcf_variant(annotvcf, source_list, source_fields):
    for line in file_util.gzopen(annotvcf):
        if annotvcf.endswith('.gz'):
            line = line.decode('UTF-8')
        if line[0] != "#":
            arr = line.split('\t')
            info = arr[7].strip()
            varkey = arr[0] + ":" + arr[1] + "_" + arr[3] + ">" + arr[4]
            var_annotfield_list = []
            print (varkey)
            for infofield in info.split(';'):
                if '=' in infofield:
                    arrf = infofield.split('=')
                    # print(arrf)
                    for section in arrf[1].split(','):
                        # print(section)
                        valuelist = section.split('|')
                        if arrf[0] in source_list:
                            if len(valuelist) != len(source_fields[arrf[0]].keys()):
                                errmsg = "\tERROR: not matching field number of "
                                errmsg += str(len(valuelist)) + " != " + str(len(source_fields[arrf[0]].keys()))
                                errmsg += " in " + arrf[0]
                                print (errmsg)
                    var_annotfield_list.append(arrf[0])

            for s1 in MUSTSOURCE:
                if s1 not in var_annotfield_list:
                    print ("\tERROR: no " + s1)
            # break


def check_annotvcf(annotvcf, mapping_table_file):
    mapping_table_col = file_util.read_table_col(mapping_table_file)
    mapping_table = file_util.read_table_key(mapping_table_file, 1)
    # print(mapping_table.keys())
    # print(mapping_table_col.keys())

    source_list, source_version_map, source_fields = get_source_list(mapping_table_col, mapping_table)
    # source_field_list = get_source_field_list(mapping_table_col, mapping_table)
    # print(source_fields['GNOMAD'].keys())
    # print(source_fields['GNOMAD']['gnomad_ac'].keys())

    # print (source_list)
    # print (source_fields['novoCaller'].keys())
    # print (source_version_map)

    check_annotvcf_header(annotvcf, source_list, source_version_map, source_fields)
    check_annotvcf_variant(annotvcf, source_list, source_fields)

    

if __name__ == "__main__":
    import proc_util
    import file_util
    MUSTSOURCE = ["VEP"]
    IGNORESOURCE = ["MUTANNO"]
    print("#USAGE: python check_annotvcf.py [ANNOT.VCF] [mapping_table.txt]")
    annotvcf = sys.argv[1]
    mapping_table_file = sys.argv[2]
    check_annotvcf(annotvcf, mapping_table_file)
