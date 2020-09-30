#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### annotviewer.oy
#### made by Daniel Minseok Kwon
#### 2020-09-29 12:39:01
#########################
import sys
import os
from ims import file_util
from ims import proc_util
import tabix

def get_first_pos(vcf):
    p = {}
    for line in file_util.gzopen(vcf):
        line = file_util.decodeb(line)
        if line[0] != "#":
            arr = line.split('\t')
            print(arr)
            p['chrom'] = arr[0]
            p['spos'] = int(arr[1])
            p['epos'] = int(arr[1])
            p['str'] = p['chrom'] + ':' + str(p['spos']) + '-' + str(p['epos'])
    return p

def pars_pos(pos_str, vcf):
    p = {}

    if pos_str == "":
        p = get_first_pos(vcf)
    else:
        arrpos = pos_str.split(':')
        p['chrom'] = arrpos[0]
        if '-' in arrpos[1]:
            arr2 = arrpos[1].split('-')
            p['spos'] = int(arr2[0])
            p['epos'] = int(arr2[1])
        else:
            p['spos'] = int(arrpos[1])
            p['epos'] = int(arrpos[1])
        p['str'] = p['chrom'] + ':' + str(p['spos']) + '-' + str(p['epos'])
    return p


def save_line(chrom, pos, ref, alt, info, fieldinfo):
    m = {}
    for field in info.split(';'):
        if "=" in field:
            arr = field.split('=')
            fieldname = arr[0]
            sections = arr[1].split(',')
            for sidx, section in enumerate(sections):
                fieldvalue = section.split('|')
                scont = ""
                if len(sections) > 1:
                    scont = " section " + str(sidx+1)

                print(fieldname + ": [" + str(len(fieldvalue)) + "] " + scont)
                
                m2 = {}
                for idx, f1 in enumerate(fieldvalue):
                    cont = '\t' + str(idx+1) + ' ' 
                    # cont += fieldname + '>' + fieldinfo[fieldname][idx] +  " : " + f1
                    # print(cont)
                    m2[fieldinfo[fieldname][idx]] = f1
                try:
                    m[fieldname]
                except KeyError:
                    m[fieldname] = []
                m[fieldname].append(m2)
                pass
        else:
            # print(field)
            pass

    # print(m['VEP'])

    cont = ""
    
    # for i in range(len(fieldinfo["VEP"])):
    #     f1 = fieldinfo["VEP"][i]
    #     cont = []
    #     cont.append(f1)
    #     for sidx in range(len(m['VEP'])):
    #         # print(m['VEP'][sidx][f1])
    #         cont.append(m['VEP'][sidx][f1])
    #     print('\t'.join(cont))

    print ("\t".join(fieldinfo["VEP"]))
    for sidx in range(len(m['VEP'])):
        cont = []
        for i in range(len(fieldinfo["VEP"])):
            f1 = fieldinfo["VEP"][i]
            cont.append(m['VEP'][sidx][f1])
        print('\t'.join(cont))
            
        
        
def get_fieldinfo_from_header(vcf):
    finfo = {}
    for line in file_util.gzopen(vcf):
        line = file_util.decodeb(line)
        if line[0] == '#':
            # print(line)
            if "##INFO=<ID=" in line:
                fieldname = line.split(',')[0].replace('##INFO=<ID=', '')
                if 'Format:' in line:
                    fieldlist = line.split('Format:\'')[-1].strip().replace('\'">', '').split('|')
                else:
                    fieldlist = ['']
                finfo[fieldname] = fieldlist
        else:
            break
    return finfo


def annotviewer(vcf, pos_str):
    pos1 = pars_pos(pos_str, vcf)

    fieldinfo = get_fieldinfo_from_header(vcf)

    tb = tabix.open(vcf)
    print(pos1['str'])

    if pos_str == "":
        for line in file_util.gzopen(vcf): 
            line = file_util.decodeb(line)
            if line[0] != '#':
                r1 = line.split('\t')
                r1[-1] = r1[-1].strip()
                chrom = r1[0]
                pos = r1[1]
                ref = r1[3]
                alt = r1[4]
                save_line(chrom, pos, ref, alt, r1[7], fieldinfo)
                break
    else:
        for r1 in tb.querys(pos1['str']): 
            chrom = r1[0]
            pos = r1[1]
            ref = r1[3]
            alt = r1[4]
            save_line(chrom, pos, ref, alt, r1[7], fieldinfo)


if __name__ == "__main__":
    print("# USAGE: python annotviewer.py annot.vcf(mti) 12:123456-123459(12:123456)")
    pos_str = ""
    if len(sys.argv) == 3:
        pos_str = sys.argv[2]
    annotviewer(sys.argv[1], pos_str)
