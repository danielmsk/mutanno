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

def load_clinvarvcf(clinvarvcf):
    global VCFCOL

    cvar = {}

    out = clinvarvcf.replace('.vcf.gz', '') + '.tsv'
    f = open(out, 'w')
    f.write('#'+'\t'.join(VCFCOL) + '\n')
    for line in file_util.gzopen(clinvarvcf):
        line = file_util.decodeb(line)
        if line[0] != '#':
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
            for f1 in VCFCOL:
                cont.append(m[f1])
            f.write('\t'.join(cont) + '\n')
            cvar[m['VARIATIONID']] = [m['CHROM'],m['POS'],m['VARIATIONID'],m['REF'],m['ALT']]
            # break
    f.close()
    print('Saved', out, len(cvar.keys()))
    return cvar

def s02_preproc_clinvarxml(clinvarxml, clinvarvcf):
    cvar = load_clinvarvcf(clinvarvcf)

    out = clinvarvcf.replace('.vcf.gz', '') + '_submission.tsv'
    f = open(out, 'w')
    cont = ['#CHROM','POS','VARIATIONID','REF','ALT', '\t'.join(XMLCOL)]
    f.write('\t'.join(cont) + '\n')
    flag = False
    i = 0
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
                    xmlfield.append(m[f1])
                cont = '\t'.join(vpos) + '\t' + '\t'.join(xmlfield)
                f.write(cont + '\n')
                # print(cont)
            except KeyError:
                pass


            flag = False
            if i % 10000 == 0:
                print(i, m)
            #     break
    
    f.close()
    print("Saved", out)
if __name__ == "__main__":
    import proc_util
    import file_util
    path = "/home/mk446/bio/mutanno/DATASOURCE/VARIANTDB/CLINVAR/hg38/"
    clinvarxml = path + "ClinVarFullRelease_2020-03.xml.gz"
    clinvarvcf = path + "clinvar_20200329.vcf.gz"
    tmpxml = "/home/mk446/bio/mutanno/SRC/scripts/clinvar/789256.xml"
    tmpxml = "/home/mk446/bio/mutanno/SRC/scripts/clinvar/568195.xml"
    s02_preproc_clinvarxml(clinvarxml, clinvarvcf)
