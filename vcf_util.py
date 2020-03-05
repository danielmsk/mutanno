#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vcf_util.py
# made by Daniel Minseok Kwon
# 2019-11-06 03:43:19
#########################


def get_info_header(infoid, infodesc, subfields, sourcesubembed=""):
    header = ""
    header += "##INFO=<ID=" + infoid + ",Number=.,Type=String"
    header += ",Description=\"" + infodesc + ". "
    if sourcesubembed != '':
        header += "Subembedded:'" + sourcesubembed + "':"
    header += "Format:'" + '|'.join(subfields) + "'\">\n"
    return header


# numgt: number format of genotype (0/0, 1/1)
# ref: reference allele
# alt: alternative allele or allels
def get_genotype(numgt, ref, alt):
    numgt = numgt.replace('|', '/')
    allele = [ref]
    allele.extend(alt.split(','))
    gt = []
    for ngt in numgt.split('/'):
        if ngt == '.':
            gt.append('.')
        else:
            gt.append(allele[int(ngt)])
    return '/'.join(gt)


def split_multiallelic_variants(vcfrecord):
    arralt = vcfrecord[4].split(',')
    rst = []
    for k in range(len(arralt)):
        r1 = []
        for j in range(len(vcfrecord)):
            if j >= 9:
                arr = vcfrecord[j].split(':')
                # GT
                phase_hipen = arr[0][1]
                if arr[0] != '0/0' and arr[0] != '0|0':
                    gt = arr[0].split('/')
                    if str(k+1) in gt:
                        if "0" in gt:
                            arr[0] = "0" + phase_hipen + "1"
                        else:
                            arr[0] = "1" + phase_hipen + "1"
                    else:
                        arr[0] = "0" + phase_hipen + "0"
                # AD
                ad = arr[1].split(',')
                arr[1] = ad[0] + ',' + ad[k+1]
                r1.append(':'.join(arr))
            else:
                r1.append(vcfrecord[j])
        r1[4] = arralt[k]

        # info
        info = []
        for f1 in r1[7].split(';'):
            if "=" in f1:
                arr = f1.split('=')
                if arr[0] in ['AC','AF','MLEAC','MLEAF']:
                    arr2 = arr[1].split(',')
                    f1 = arr[0] + "=" + arr2[k]
            info.append(f1)
        info.append('multiallele=' + vcfrecord[0] + ':' + vcfrecord[1] + ' ' + vcfrecord[3] + '/' + vcfrecord[4])
        r1[7] = ';'.join(info)
        rst.append(r1)
    return rst


def encode_value(v1):
    v1 = v1.replace(';', "%3B").replace('=', "%3D").replace("|", "%7C").replace(",", "%2C")
    v1 = v1.replace('"', "%22").replace('~', "%7E").replace(' ', '%20')
    return v1


def encode_infovalue(v1, delimiter=""):
    if v1 == '.' or v1 == '-':
        v1 = ''
    v1 = encode_value(v1)
    if delimiter != "":
        v1 = v1.replace(encode_value(delimiter), '~')
    # v1 = urllib.parse.quote(v1)
    return v1


def remove_nonvariant_base(ref, alt, pos=0, miss_char=''):
    p = 0
    if min(len(ref), len(alt)) >= 1:
        for j in range(min(len(ref), len(alt))):
            if len(ref) == 0 or len(alt) == 0:
                break
            if alt[0] == ref[0]:
                alt = alt[1:]
                ref = ref[1:]
                p += 1
            else:
                break
        for i in range(1, min(len(ref), len(alt))):
            if len(ref) == 0 or len(alt) == 0:
                break
            if alt[-1] == ref[-1]:
                alt = alt[:-1]
                ref = ref[:-1]
            else:
                break
    if ref == '':
        ref = miss_char
    if alt == '':
        alt = miss_char
    pos += p
    return (ref, alt, pos)


class VCFHEADER():
    collist = []
    mutannofield = {}
    header = ""

    def __init__(self, header):
        self.header = header
        self.pars()

    def pars(self):
        flag_mutanno = False
        m = {}
        for line in self.header.split('\n'):
            if line[:2] == "##":
                if line[:len("##MUTANNO=<ID=")] == "##MUTANNO=<ID=":
                    arr = line.replace('##MUTANNO=<ID=', '').replace('>', '').split(',')
                    d = {}
                    d['sourcename'] = arr[0].strip()
                    for f1 in arr[1:]:
                        arr2 = f1.strip().replace('""', '').split('=')
                        d[arr2[0].strip()] = arr2[1].strip()
                    d['fieldlist'] = []
                    m[arr[0].strip()] = d

                if line[:len("##MUTANNO=<ID=MUTANNO")] == "##MUTANNO=<ID=MUTANNO":
                    flag_mutanno = True

                if line[:len("##INFO=")] == "##INFO=" and 'Format:' in line and flag_mutanno:
                    arr = line.strip().split(',')
                    sourcename = arr[0].replace('##INFO=<ID=', '')
                    if sourcename in m.keys():
                        m[sourcename]['fieldlist'] = line.split(
                            'Format:')[-1].replace("'", "").replace('">', '').strip().split('|')
            elif line[:len("#CHROM")] == "#CHROM":
                self.collist = line[1:].strip().split('\t')
        self.mutannofield = m


#     fieldkeynamelist.append( conv_cgap_filed(k1 + '.' + v1))
# f.write('\t'.join(fieldkeynamelist)+'\n')

def pars_vcfline(vcfline, collist=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAT', 'FILTER', 'INFO', 'FORMAT']):
    arr = vcfline.split('\t')
    arr[-1] = arr[-1].strip()

    sid = 0
    v1 = {}
    for i in range(len(arr)):
        try:
            colname = collist[i]
        except IndexError:
            sid += 1
            colname = 'S' + str(sid)

        if colname == "FORMAT":
            sformat = arr[i].split(':')

        if colname == "INFO":
            d = {}
            for f1 in arr[i].split(';'):
                f1 = f1.strip()
                if '=' in f1:
                    arr2 = f1.split('=')
                    d[arr2[0]] = arr2[1]
                else:
                    d[f1] = f1
        if i >= 9:
            d = {}
            arr2 = arr[i].split(':')
            for j in range(len(arr2)):
                d[sformat[j]] = arr2[j]

        if i >= 9 or colname == "INFO":
            v1[colname] = d
        else:
            v1[colname] = arr[i]
    return v1
