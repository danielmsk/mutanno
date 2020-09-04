#!/usr/bin/env python
# -*- coding: utf-8 -*-

from . import struct_util
import math
from itertools import combinations_with_replacement
VCF_COL = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLESTART']

INFOIDX = VCF_COL.index('INFO')

# https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39
CHROM_HGVS = {}
CHROM_HGVS["1"] = "NC_000001.11"
CHROM_HGVS["2"] = "NC_000002.12"
CHROM_HGVS["3"] = "NC_000003.12"
CHROM_HGVS["4"] = "NC_000004.12"
CHROM_HGVS["5"] = "NC_000005.10"
CHROM_HGVS["6"] = "NC_000006.12"
CHROM_HGVS["7"] = "NC_000007.14"
CHROM_HGVS["8"] = "NC_000008.11"
CHROM_HGVS["9"] = "NC_000009.12"
CHROM_HGVS["10"] = "NC_000010.11"
CHROM_HGVS["11"] = "NC_000011.10"
CHROM_HGVS["12"] = "NC_000012.12"
CHROM_HGVS["13"] = "NC_000013.11"
CHROM_HGVS["14"] = "NC_000014.9"
CHROM_HGVS["15"] = "NC_000015.10"
CHROM_HGVS["16"] = "NC_000016.10"
CHROM_HGVS["17"] = "NC_000017.11"
CHROM_HGVS["18"] = "NC_000018.10"
CHROM_HGVS["19"] = "NC_000019.10"
CHROM_HGVS["20"] = "NC_000020.11"
CHROM_HGVS["21"] = "NC_000021.9"
CHROM_HGVS["22"] = "NC_000022.11"
CHROM_HGVS["X"] = "NC_000023.11"
CHROM_HGVS["Y"] = "NC_000024.10"


def parse_info_header_line(line, d={}):
    fields = line.replace('##INFO=<', '').split(',')
    d2 = {}
    for idx, f1 in enumerate(fields):
        arr = f1.split('=')
        if arr[0] == "Description":
            if idx+1 == len(fields):
                desc = arr[1]
            else:
                desc = arr[1] + ',' + ",".join(fields[(idx+1):])
            desc = desc.strip()[:-1].replace('"', '')
            d2[arr[0]] = desc
            if 'Format:' in desc:
                arr2 = desc.split('Format:')
                d2['Format'] = arr2[-1].replace("'", "").split('|')
                if 'Subembedded:' in desc:
                    arr3 = arr2[0].split('Subembedded:')
                    d2['Subembedded'] = arr3[-1].strip()
            break
        else:
            d2[arr[0]] = arr[1]
    d[d2['ID']] = d2
    return d


def get_hgvsg(chrom, pos, ref, alt):
    try:
        hgvsg = CHROM_HGVS[chrom.replace('chr', '')] + ':g.'
        if len(ref) > len(alt):
            if len(ref) == 2:
                hgvsg += str(pos + 1) + 'del'
            else:
                hgvsg += str(pos + 1) + '_' + str(pos + len(ref) - 1) + 'del'
        elif len(ref) < len(alt):
            hgvsg += str(pos) + '_' + str(pos + 1) + 'ins' + alt[1:]
        else:
            hgvsg += str(pos) + ref + '>' + alt
    except KeyError:
        hgvsg = ""
    return hgvsg


def get_variant_class(ref, alt):
    # TODO: need to upgrade for ',' in alt.
    vcls = ""
    if len(ref) == len(alt) and len(alt) == 1:
        vcls = "SNV"
    elif len(ref) < len(alt):
        vcls = "INS"
    elif len(ref) > len(alt):
        vcls = "DEL"
    return vcls


def add_info(info1, info2):
    info1 = strip_info(info1)
    info2 = strip_info(info2)
    if info1 == "":
        info = info2
    if info2 == "":
        info = info1
    if info1 != "" and info2 != "":
        info = info1 + ';' + info2
    return info


def strip_info(info):
    if info == ".":
        info = ""
    if info != "" and info[-1] == ";":
        info = info[:-1]
    return info


def convert_to_metadata(d):
    cont = []
    for field in d.keys():
        if isinstance(d[field], dict):
            for key in d[field].keys():
                valuearr = ["ID=" + key]
                for k2 in d[field][key].keys():
                    if k2 in ['Number', 'Type']:
                        valuearr.append(k2 + '=' + d[field][key][k2])
                    else:
                        valuearr.append(k2 + '="' + d[field][key][k2] + '"')
                value = '<' + ','.join(valuearr) + '>'
                cont.append('##'+field + '=' + value)
        else:
            value = d[field]
            cont.append('##'+field + '=' + value)
    return '\n'.join(cont)


def get_info_header(headertype, infoid, version, vdate,  infodesc, subfields, sourcesubembed=""):
    header = ""

    if headertype == "MUTANNO":
        header += "##MUTANNO=<ID="+infoid+",Version=\""+version+"\",Date=\""+vdate+"\">\n"
    elif headertype == "INFO":
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


def get_numgt(gt, ref, alt, delimiter=''):
    if delimiter == '':
        if '|' in gt:
            delimiter = '|'
        else:
            delimiter = '/'
        gt = gt.replace(delimiter, '/')

    numgt = []
    for g1 in gt.split('/'):
        if g1 == alt:
            numgt.append('1')
        else:
            numgt.append('0')

    return delimiter.join(numgt)


def split_multiallelic_variants(vcfrecord):
    # print('>vcfrecord:', vcfrecord)
    arralt = vcfrecord[4].split(',')
    arrformat = vcfrecord[8].split(':')
    rst = []
    for k in range(len(arralt)):
        r1 = []
        ac = 0
        for j in range(len(vcfrecord)):
            if j >= 9:
                arr = vcfrecord[j].strip().split(':')
                # GT
                phase_hipen = arr[0][1]
                if arr[0] == "."+phase_hipen+".":
                    pass
                elif arr[0] != '0' + phase_hipen + '0':
                    gt = arr[0].split(phase_hipen)
                    if str(k+1) in gt:
                        if "0" in gt:
                            arr[0] = "0" + phase_hipen + "1"
                            ac += 1
                        else:
                            arr[0] = "1" + phase_hipen + "1"
                            ac += 2
                    else:
                        arr[0] = "0" + phase_hipen + "0"
                # AD
                ad = arr[arrformat.index('AD')].split(',')
                arr[arrformat.index('AD')] = ad[0] + ',' + ad[k+1]

                # DP
                if arrformat.index('DP') < len(arr):
                    arr[arrformat.index('DP')] = str(int(ad[0]) + int(ad[k+1]))

                # PL
                if arrformat.index('PL') < len(arr):
                    arr_pl = arr[arrformat.index('PL')].strip().split(',')
                    if len(arr_pl) > 1:
                        new_pl = get_biallelepl_multiallelepl('0/' + str(k+1), arr_pl)
                        arr[arrformat.index('PL')] = ','.join(new_pl)
                        
                r1.append(':'.join(arr))
            else:
                r1.append(vcfrecord[j])
        r1[4] = arralt[k]

        # info
        info = []
        for f1 in r1[7].split(';'):
            if "=" in f1:
                arr = f1.split('=')
                if arr[0] == 'AC':
                    f1 = arr[0] + "=" + str(ac)
                
                if arr[0] == 'AF':
                    sample_size = (len(vcfrecord) - 9)
                    f1 = arr[0] + "=" + str(round(ac / (sample_size * 2), 3))

                ### if we want to use AC, AF from VCF,
                # if arr[0] in ['AC', 'AF', 'MLEAC', 'MLEAF']:
                if arr[0] in ['MLEAC', 'MLEAF']:
                    arr2 = arr[1].split(',')
                    f1 = arr[0] + "=" + arr2[k]

            info.append(f1)
        # info.append('multiallele=' + vcfrecord[0] + ':' + vcfrecord[1] + '%20' + vcfrecord[3] + '/' + vcfrecord[4])
        r1[7] = ';'.join(info)
        rst.append(r1)

    return rst


def get_biallelepl_multiallelepl(numgt, multiallelepl):
    if '|' in numgt:
        delimiter = '|'
    else:
        delimiter = '/'

    # from F(j/k) = (k*(k+1)/2)+j of https://samtools.github.io/hts-specs/VCFv4.1.pdf
    no_allele = int( (math.sqrt(1 + 8 * len(multiallelepl)) - 1) / 2)
    comb_gt = list(combinations_with_replacement(list(range(no_allele)), 2))

    pl_pos_map = {}
    for cgt in comb_gt:
        gt = '/'.join(struct_util.convto_str_array(cgt))
        # from F(j/k) = (k*(k+1)/2)+j of https://samtools.github.io/hts-specs/VCFv4.1.pdf
        pl_pos_map[gt] = int((cgt[1] * (cgt[1] + 1) / 2) + cgt[0])

    new_pl = []
    for ngt in list(combinations_with_replacement(numgt.split(delimiter), 2)):
        try:
            pidx = pl_pos_map['/'.join(ngt)]
            new_pl.append(multiallelepl[pidx])
        except KeyError:
            print('KeyError:',numgt, multiallelepl, pl_pos_map, ngt)
            new_pl.append('0')
    return new_pl


def get_numgt_from_multiallele(genotype, multigenotype):
    if '|' in genotype:
        delimiter = '|'
    else:
        delimiter = '/'
    mgt1 = multigenotype.split('/')
    m = {}
    m[mgt1[0]] = '0'
    for i, alt in enumerate(mgt1[1].split(',')):
        m[alt] = str(i+1)
    numgt = []
    for b1 in genotype.split(delimiter):
        numgt.append(m[b1])
    return delimiter.join(numgt)


def encode_value(v1):
    v1 = v1.replace(';', "%3B").replace('=', "%3D").replace("|", "%7C").replace(",", "%2C")
    v1 = v1.replace('"', "%22").replace('~', "%7E").replace(' ', '%20')
    return v1


def decode_value(v1):
    v1 = v1.replace("%3B", ';').replace("%3D", '=').replace("%7C", "|").replace("%2C", ",")
    v1 = v1.replace("%22", '"').replace("%7E", '~').replace('%20', ' ')
    return v1


def encode_infovalue(v1, delimiter=""):
    if v1 is None:
        rst = ""
    elif isinstance(v1, bool):
        if v1:
            rst = "1"
        else:
            rst = "0"
    else:
        rst = str(v1)
        if rst == '.' or rst == '-':
            rst = ''
        rst = encode_value(rst)
        if delimiter != "":
            rst = rst.replace(encode_value(delimiter), '~')

    return rst


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


def pars_info_header(infoheader):
    d = {}
    for s1 in infoheader.strip().split(";"):
        if "=" in s1:
            arr = s1.split('=')
            d[arr[0]] = arr[1].strip().split('|')
    return(d)


def pars_info_field(infofield):
    d = {}
    for s1 in infofield.strip().split(";"):
        if "=" in s1:
            arr = s1.split('=')
            attr = []
            for sec in arr[1].split(','):
                attr.append(sec.strip().split('|'))
            d[arr[0]] = attr
        else:
            d[s1] = True
    return(d)


def pars_info_field_with_infoheader(infofield, infoheaderdict):
    d = {}
    for s1 in infofield.strip().split(";"):
        if "=" in s1:
            arr = s1.split('=')
            fname = arr[0]
            attr = []
            for sec in arr[1].split(','):
                d2 = {}
                fields = sec.strip().split('|')

                for idx, f1 in enumerate(infoheaderdict[fname]):
                    d2[f1] = fields[idx]

                attr.append(d2)
            d[fname] = attr
        else:
            d[s1] = True
    return(d)

