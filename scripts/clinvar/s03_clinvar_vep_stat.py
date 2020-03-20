#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s12_clinvar_vep_stat.py
# made by Daniel Minseok Kwon
# 2020-02-10 15:29:47
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
import tabix

class MutAnnoDataTSV():
    def __init__(self, filepath):
        self.file = filepath
        self.tbchrom = ""
        self.tb = ""
        self.colidx = {}
        self.set_header()

    def load_tabixpointer(self):
        self.tb = tabix.open(self.file.replace('#CHROM#',self.tbchrom))

    def get_annot_col(self, colcont):
        m = {}
        for f1 in colcont.split(';'):
            arr1 = f1.split('=')
            m[arr1[0]] = arr1[1].split('|')
        return m


    def set_header(self):
        print(self.file.replace('#CHROM#','1'))
        for line in file_util.gzopen(self.file.replace('#CHROM#','1')):
            line = file_util.decodeb(line)
            if line[:len('#CHROM')] == '#CHROM':
                arr = line[1:].split('\t')
                for j in range(len(arr[:-1])):
                    self.colidx[arr[j]] = j
                self.annotcol = self.get_annot_col(arr[-1].strip())
                break

    def get_variant(self, chrom, pos, ref="", alt=""):
        if self.tbchrom != chrom:
            self.tbchrom = chrom
            self.load_tabixpointer()

        mvar_list = []
        for rec in self.tb.querys(chrom + ":"+str(pos)+"-"+str(pos)):
            if rec[self.colidx['CHROM']] == chrom and rec[self.colidx['POS']] == str(pos):
                if (ref == "" or rec[self.colidx['REF']] == ref) and (alt == "" or rec[self.colidx['ALT']] == alt):
                    mvar_list.append(MutAnnoDataVariant(rec, self.annotcol, self.colidx))
        return mvar_list

class MutAnnoDataVariant():
    def __init__(self, record, fieldnames, colidx):
        self.fieldnames = fieldnames
        self.colidx = colidx
        self.annot = {}
        self.record = record
        self.pars_record(record)
    
    def pars_record(self, record):
        self.chrom = record[self.colidx['CHROM']]
        self.pos = int(record[self.colidx['POS']])
        self.ID = record[self.colidx['ID']]
        self.ref = record[self.colidx['REF']]
        self.alt = record[self.colidx['ALT']]

        for a1 in record[-1].split(';'):
            a2 = a1.split('=')
            skey = a2[0]
            self.annot[skey] = []
            for a3 in a2[1].split(','):
                a4 = a3.split('|')
                d = {}
                for i in range(len(self.fieldnames[skey])):
                    d[self.fieldnames[skey][i]] = a4[i]
                self.annot[skey].append(d)

class CountStat():
    def __init__(self):
        self.cnt = {}

    def add(self, key1, key2, key3=''):
        if key1 not in self.cnt.keys():
            self.cnt[key1] = {}

        if key3 == '':
            try:
                self.cnt[key1][key2] += 1
            except KeyError: 
                self.cnt[key1][key2] = 1
        else:
            if key2 not in self.cnt[key1].keys():
                self.cnt[key1][key2] = {}

            try:
                self.cnt[key1][key2][key3] += 1
            except KeyError: 
                self.cnt[key1][key2][key3] = 1




def s12_clinvar_vep_stat(clinvartsv, vepsource, out):

    mannot = MutAnnoDataTSV(mutannot)

    cstat = CountStat()
    f = open(out, 'w')
    i = 0
    for line in file_util.gzopen(clinvartsv):
        line = line.decode('UTF-8')
        if line[0] != '#':
            i += 1
            arr = line.split('\t')
            # print(arr)
            chrom = arr[0]
            pos = int(arr[1])
            ref = arr[2]
            alt = arr[3]
            clnsig = arr[11]
            cstat.add('CLNSIG', clnsig)

            varannotlist = mannot.get_variant(chrom, pos, ref, alt)
            if len(varannotlist) == 1:
                varannot = varannotlist[0]
                # print(varannot.fieldnames)
                if 'VEP' in varannot.annot.keys():
                    # print(chrom, pos, ref, alt)
                    print(chrom, pos, ref, alt, len(varannot.annot['VEP']))
                    for vep in varannot.annot['VEP']:
                        # print(vep['Consequence'])
                        cstat.add('CLNSIG2', clnsig, vep['Consequence'])
                    # print(clnsig, vep['Consequence'], vep['Gene'])
                # print(clnsig)
            if i > 1000:
                break

    for k1 in cstat.cnt.keys():
        for k2 in cstat.cnt[k1].keys():
            print(k1, k2, cstat.cnt[k1][k2])
    f.close()

if __name__ == "__main__":
    import proc_util
    import file_util
    import vcf_util
    clinvartsv = "/home/mk446/mutanno/DATASOURCE/VARIANTDB/CLINVAR/hg38/clinvar.20200106.hg38.tsv.gz"
    mutannot = "/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/v0.3.2_200309/mvp_datasource_v0.3.2_200309.chr#CHROM#.tsi.gz"
    out = clinvartsv.replace(".tsv.gz",".tsv") + ".vep.stat"
    s12_clinvar_vep_stat(clinvartsv, mutannot, out)
