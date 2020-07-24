#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s02_preproc_clinvarxml.py
# made by Daniel Minseok Kwon
# 2020-05-05 09:47:16
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
# from xml.dom.minidom import parse, parseString
import xmltodict
import json

VCFCOL = ["CHROM","POS","VARIATIONID","REF","ALT","ALLELEID","CLNDN","CLNDNINCL","CLNDISDB","CLNDISDBINCL","CLNHGVS","CLNREVSTAT","CLNSIG","CLNSIGCONF","CLNSIGINCL","CLNVC","CLNVCSO","CLNVI","GENEINFO","MC","ORIGIN","SSR"]
# XMLCOL = ["ClinVarAccession","Interpretation","DateLastEvaluated","ReviewStatus","Method","Condition","AlleleOrigin","Submitter","SubmitterID","Citation","Comment"]
XMLCOL = ["ClinVarAccession","Interpretation","DateLastEvaluated","ReviewStatus","Method","Condition","AlleleOrigin","Submitter","SubmitterID","Citation"]


def get_dict_data(dict1, keylist):
    rst = ""
    rst1 = dict1
    for i, key1 in enumerate(keylist):
        # try:
        #     rst1[key1]
        # except TypeError:
        #     print('TypeError:',key1, rst1)
        if type(rst1) == list:
            rst1 = rst1[0]

        try:
            rst1[key1]
            rst1 = rst1[key1]
            rst = rst1
        except KeyError:
            rst = ""
            break
    return rst

def cnt_field(cnt, f1, v1):
    try:
        cnt[f1]
    except KeyError:
        cnt[f1]={}
    try:
        cnt[f1][v1] += 1
    except KeyError:
        cnt[f1][v1] = 1
    return cnt

CONV_SIGN_LIST = []
CONV_SIGN_LIST.append((';', "%3B"))
CONV_SIGN_LIST.append(('=', "%3D"))
CONV_SIGN_LIST.append(("|", "%7C"))
CONV_SIGN_LIST.append((",", "%2C"))
CONV_SIGN_LIST.append(('"', "%22"))
CONV_SIGN_LIST.append(('~', "%7E"))
CONV_SIGN_LIST.append((' ', '%20'))

def encode_value(v1, except_list=[]):
    for s1 in CONV_SIGN_LIST:
        if s1[0] not in except_list:
            v1 = v1.replace(s1[0], s1[1])
    # v1 = v1.replace(';', "%3B").replace('=', "%3D").replace("|", "%7C").replace(",", "%2C")
    # v1 = v1.replace('"', "%22").replace('~', "%7E").replace(' ', '%20')
    return v1

# 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other


def load_clinvarvcf(clinvarvcf):
    global VCFCOL

    i = 0
    cvar = {}
    cnt = {}
    out = clinvarvcf.replace('.vcf.gz', '') + '.tsi'
    f = open(out, 'w')
    f.write('#'+'\t'.join(VCFCOL[:5]) + '\tCLINVAR=' + '|'.join(VCFCOL[5:]) + '\n')
    for line in file_util.gzopen(clinvarvcf):
        line = file_util.decodeb(line)
        if line[0] != '#':
            i += 1
            arr = line.split('\t')
            m = {}
            for f1 in VCFCOL:
                m[f1] = ""
            m['CHROM'] = arr[0].strip()
            if m['CHROM'] == "MT":
                m['CHROM'] = "M"
            m['POS'] = arr[1].strip()
            m['VARIATIONID'] = arr[2].strip()
            m['REF'] = arr[3].strip()
            m['ALT'] = arr[4].strip()
            for f1 in arr[7].strip().split(';'):
                arr2 = f1.split('=')
                m[arr2[0]] = arr2[1].strip()
            
            cont = []
            for f1 in VCFCOL[:5]:
                cont.append(m[f1])

            info = []
            for idx in range(5, len(VCFCOL)):
                f1 = VCFCOL[idx]
                if f1 == "MC":
                    arr2 = m[f1].split(',')
                    arr3 = []
                    for a2 in arr2:
                        arr3.append(a2.split('|')[0])
                    m[f1] = '~'.join(arr3)
                elif f1 == "CLNDISDB":
                    m[f1] = m[f1].replace('.|','').replace('.|','').replace('.|','').replace('.|','').replace('.|','').replace('.|','')
                    m[f1] = m[f1].replace('.|','').replace('.|','').replace('.|','').replace('.|','').replace('.|','').replace('.|','')
                    m[f1] = m[f1].replace('.','')
                    m[f1] = m[f1].replace('|','~').replace(',','~')
                    m[f1] = encode_value(m[f1], ['~'])
                elif f1 == "CLNDN":
                    m[f1] = m[f1].replace('|','~')
                    m[f1] = encode_value(m[f1], ['~'])
                elif f1 == "CLNDN":
                    m[f1] = m[f1].replace('|','~')
                    m[f1] = encode_value(m[f1], ['~'])
                elif f1 == "CLNSIGINCL":
                    m[f1] = m[f1].replace('|','~')
                    m[f1] = encode_value(m[f1], ['~'])
                elif f1 == "CLNSIGCONF":
                    m[f1] = m[f1].replace(',','~')
                    m[f1] = encode_value(m[f1], ['~'])
                elif f1 == "CLNVI":
                    m[f1] = m[f1].replace('|','~')
                    m[f1] = encode_value(m[f1], ['~'])
                elif f1 == "CLNDISDBINCL":
                    m[f1] = m[f1].replace('.|','').replace('.|','').replace('.|','').replace('.|','').replace('.|','').replace('.|','')
                    m[f1] = m[f1].replace('.|','').replace('.|','').replace('.|','').replace('.|','').replace('.|','').replace('.|','')
                    m[f1] = m[f1].replace('.','')
                    m[f1] = m[f1].replace('|','~').replace(',','~')
                    m[f1] = encode_value(m[f1], ['~'])
                elif f1 == "GENEINFO":
                    m[f1] = m[f1].replace('|','~')
                    m[f1] = encode_value(m[f1], ['~'])
                else:
                    m[f1] = encode_value(m[f1])
                info.append(m[f1])

                cnt = cnt_field(cnt, f1, m[f1])

                
            cont.append('CLINVAR='+'|'.join(info))
            f.write('\t'.join(cont) + '\n')
            cvar[m['VARIATIONID']] = [m['CHROM'],m['POS'],m['VARIATIONID'],m['REF'],m['ALT']]
            # break

            if i % 100000 == 0:
                print(i, m)
                # break

    f.close()
    print('Saved', out, len(cvar.keys()))

    for k1 in VCFCOL[5:]:
        statfile = out + "_stat/stat_" + k1.replace(' ','_') + '.cnt'
        file_util.check_dir(statfile)
        cont = ''
        ks = cnt[k1].keys()
        for k2 in sorted(ks):
            # print(k1, k2, cnt[k1][k2])
            cont += k2 + '\t' + str(cnt[k1][k2]) + '\n'
        file_util.fileSave(statfile, cont, 'w')

    return cvar

def s02_preproc_clinvarxml(clinvarxml, clinvarvcf):
    cvar = load_clinvarvcf(clinvarvcf)

    out = clinvarvcf.replace('.vcf.gz', '') + '_submission.tsi'
    f = open(out, 'w')
    cont = ['#CHROM','POS','VARIATIONID','REF','ALT', "CLINVAR_SUBMISSION="+'|'.join(XMLCOL)]
    f.write('\t'.join(cont) + '\n')
    flag = False
    i = 0
    cnt = {}
    for line in file_util.gzopen(clinvarxml):
        line = file_util.decodeb(line)
        
        if '<ClinVarSet ' in line:
            flag = True
            setcont = ""
            
            
        if flag:
            setcont += line

        if '</ClinVarSet>' in line:
            i += 1
            
            xml = xmltodict.parse(setcont)
            cset = xml['ClinVarSet']
            aset = get_dict_data(xml, ['ClinVarSet', 'ClinVarAssertion'])

            # print(setcont)

            m = {}
            for f1 in XMLCOL:
                m[f1] = ''
            m['VARIATIONID'] = get_dict_data(xml,['ClinVarSet','ReferenceClinVarAssertion','MeasureSet','@ID'])
            m['ClinVarAccession'] = get_dict_data(aset, ['ClinVarAccession','@Acc'])
            m['Interpretation'] = get_dict_data(aset, ['ClinicalSignificance','Description'])
            m['DateLastEvaluated'] = get_dict_data(aset, ['ClinicalSignificance','@DateLastEvaluated'])
            m['ReviewStatus'] = get_dict_data(aset, ['ClinicalSignificance','ReviewStatus'])
            # m['Method'] = get_dict_data(aset, ['ObservedIn','Method', 'MethodType'])
            m['Method'] = get_dict_data(xml, ['ClinVarSet','ReferenceClinVarAssertion','ObservedIn','Method','MethodType'])
            m['Condition'] = get_dict_data(aset, ['TraitSet','Trait','Name','ElementValue','#text'])
            m['AlleleOrigin'] = get_dict_data(xml, ['ClinVarSet','ReferenceClinVarAssertion','ObservedIn','Sample','Origin'])
            # m['AlleleOrigin'] = get_dict_data(aset, ['ObservedIn','Sample','Origin'])
            m['Submitter'] = get_dict_data(aset, ['ClinVarSubmissionID','@submitter'])
            m['SubmitterID'] = get_dict_data(aset, ['ClinVarAccession','@OrgID']) # https://www.ncbi.nlm.nih.gov/clinvar/submitters/<ID>/
            m['Citation'] = get_dict_data(aset, ['Citation','ID','#text'])
            # m['Comment'] = get_dict_data(aset, ['ClinicalSignificance','Comment','#text'])
            # print(m)
            try:
                vpos = cvar[m['VARIATIONID']]
                xmlfield = []
                for f1 in XMLCOL:
                    m[f1] = encode_value(m[f1])
                    xmlfield.append(m[f1])
                    cnt = cnt_field(cnt, f1, m[f1])

                cont = '\t'.join(vpos) + '\t' + "CLINVAR_SUBMISSION="+ '|'.join(xmlfield)
                f.write(cont + '\n')
                # print(cont)
            except KeyError:
                pass


            flag = False
            if i % 10000 == 0:
                print(i, m)
                # break
    
    f.close()

    for k1 in XMLCOL:
        statfile = out + "_stat/stat_" + k1.replace(' ','_') + '.cnt'
        file_util.check_dir(statfile)
        cont = ''
        ks = cnt[k1].keys()
        for k2 in sorted(ks):
            # print(k1, k2, cnt[k1][k2])
            cont += k2 + '\t' + str(cnt[k1][k2]) + '\n'
        file_util.fileSave(statfile, cont, 'w')


    print("Saved", out)
if __name__ == "__main__":
    import proc_util
    import file_util
    path = "/home/mk446/bio/mutanno/DATASOURCE/VARIANTDB/CLINVAR/hg38/"
    clinvarxml = path + "ClinVarFullRelease_2020-03.xml.gz"
    clinvarvcf = path + "clinvar_20200329.vcf.gz"
    # tmpxml = "/home/mk446/bio/mutanno/SRC/scripts/clinvar/789256.xml"
    # tmpxml = "/home/mk446/bio/mutanno/SRC/scripts/clinvar/568195.xml"
    s02_preproc_clinvarxml(clinvarxml, clinvarvcf)
    # cvar = load_clinvarvcf(clinvarvcf)
