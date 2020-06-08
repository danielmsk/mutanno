
import re
import vcf
import argparse
import json
from .util import file_util
from .util import vcf_util
from . import _options
from .util.struct_util import get_dict_value as dv
from .model.datasource import DataSourceList
from . import external_functions

INFOIDX = vcf_util.VCF_COL.index('INFO')
FORMATIDX = vcf_util.VCF_COL.index('FORMAT')


def pars_header(line):
    d = {}
    p1 = re.compile(r'''##(?P<field>.+?)=(?P<val>.+)''')
    p2 = re.compile(r'''(?P<key>.+?)=(?P<val>[^,]+)(,|)''')
    m1 = p1.match(line)
    field = m1.group('field')
    value = m1.group('val')
    if value[0] != "<" and value[-1] != ">":
        d[field] = value
    else:
        for m2 in p2.finditer(value[1:-1]):
            key = m2.group('key')
            val = m2.group('val')
            try:
                d[field][key] = val
            except KeyError:
                d[field] = {}
                d[field][key] = val
    return d

def read_metadata_and_check_duplicate(header):
    meta = {}
    for line in header.strip().split('\n'):
        line = line.strip()
        d = pars_header(line)
        field = list(d.keys())[0]
        
        try:
            meta[field]
        except KeyError:
            if type(d[field]) == type({}):
                meta[field] = {}
            else:
                meta[field] = []

        if type(d[field]) == type({}):
            fid = d[field]['ID']
            try:
                meta[field][fid]
                raise Exception("VCF header duplicates: " + line)
            except KeyError:
                meta[field][fid] = d[field]
        else:
            if d[field] in meta[field]:
                raise Exception("VCF header duplicates: " + line)
            meta[field].append(d[field])
    return meta


class AnnotVCFValidator:
    def __init__(self, vcf_file):
        self.vcf_file = vcf_file
        self.reader = vcf.Reader(open(self.vcf_file))
        self.header = ""
        self.colnames = []
        self.mutanno_sources = {}
        self.dslist = None
        self.opt = None
        self.varkey = ""
        self.read_header()

    def set_datastructure(self, jsonfile=""):
        self.dslist = DataSourceList(jsonfile)

    def set_option(self, opt):
        if isinstance(opt, list):
            self.opt = _options.get_options()
        else:
            self.opt = opt

    def read_header(self):
        for line in file_util.gzopen(self.vcf_file):
            line = file_util.decodeb(line)
            if line[:2] == '##':
                self.header += line
            elif line[0] == '#':
                self.colnames = line.strip().split('\t')
            else:
                break

    def load_mutanno_structure(self):
        for line in self.reader._header_lines:
            if line[:len('##MUTANNO=<ID=')] == '##MUTANNO=<ID=':
                d = {}
                for f1 in line.replace('##MUTANNO=<','').replace('>','').split(','):
                    arr2 = f1.split('=')
                    d[arr2[0]] = arr2[1].replace('"','')
                self.mutanno_sources[d['ID']] = d

        for s1 in self.mutanno_sources.keys():
            try:
                self.reader.infos[s1]
            except KeyError:
                self.raise_exception("MUTANNO's " + infokey + " source doesn't have INFO field in header.")

        for infokey in self.reader.infos.keys():
            desc = self.reader.infos[infokey].desc
            if 'Format:' in desc:
                fields = desc.split("Format:'")[-1].replace("'","").split('|')
                try:
                    self.mutanno_sources[infokey]
                    self.mutanno_sources[infokey]['fields'] = fields
                except KeyError:
                    self.raise_exception("MUTANNO's " + infokey + " source doesn't have MUTANNO's version field in header.")

    def raise_exception(self, msg, print_varkey=True):
        if print_varkey:
            msg += " " + self.varkey
        msg += " in " + self.vcf_file
        raise Exception(msg)   

    def warn(self, msg, print_varkey=True):
        if print_varkey:
            msg += " " + self.varkey
        print("\tWARN: " + msg + " " + self.vcf_file)

    def check_vep_consequence_and_most_severe(self, info, col, v1=None):
        considx = col.index('Consequence')

        if 'VEP' in info.keys():
            vep_sections = []
            # print(info['VEP'])
            most_severe = 0
            for idx, attr in enumerate(info['VEP']):
                if isinstance(attr, str):
                    fields = attr.split('|')
                else:
                    fields = attr
                vep_sections.append(fields)
                if '&' in fields[considx]:
                    self.raise_exception("VEP Consequence has '&'. " + fields[considx])

                # print('col:',col)
                if 'MOST_SEVERE' in col:
                    most_severe += int(fields[col.index('MOST_SEVERE')])

            if 'MOST_SEVERE' in col and most_severe != 1:
                self.raise_exception("\tThis variant doesn't have MOST_SEVERE transcript. ")
            # for idx, fields in enumerate(vep_sections):
            #     is_most_severe = external_functions.is_most_severe_transcript(vep_sections, idx, self.mutanno_sources['VEP']['fields'])
            #     if str(is_most_severe) != str(fields[mostidx]):
            #         pass
            #         # raise Exception("VEP MOST_SEVERE is not matched." + v1.CHROM + ':' + str(v1.POS) + " " + fields[1] + " " + str(is_most_severe) + ":" + fields[mostidx])
    
    def check_genes(self, info):
        try:
            for f1 in info["GENES"]:
                if "tmpOOOOOOOOOO" in f1:
                    # fields = f1.split('|')
                    self.raise_exception("\tThe add_severe_consequence function of GENES fields doesn't work.")
        except KeyError:
            if 'VEP' in info.keys() and 'GENES' in self.dslist.sources.keys():
                self.raise_exception("\tNo GENES.")
            else:
                self.warn("No GENES.")

    def check_feature_ncbi(self, info):
        try:
            for vep in info["VEP"]:
                vepfields = vep.split('|')
                if "Feature_ncbi" in self.mutanno_sources['VEP']['fields']:
                    idx = self.mutanno_sources['VEP']['fields'].index("Feature_ncbi")
                    if vepfields[idx] == "":
                        self.warn("No Feature_ncbi")
                    else:
                        print("\tOK. Feature_ncbi=" + vepfields[idx])
        except KeyError:
            # self.raise_exception("\tNo VEP.")
            self.warn("No VEP.")


    def check_hg19(self, info):
        try:
            info["HG19"]
            for hg19 in info["HG19"]:
                if len(hg19.split('|')) != len(self.mutanno_sources['HG19']['fields']):
                    self.raise_exception("HG19 is unmatched " + hg19)
                else:
                    print('\tOK. HG19='+','.join(info["HG19"]))
        except KeyError:
            self.warn("\tNo HG19.")

    def check_samplegeno(self, info):
        try:
            # print(info["SAMPLEGENO"])
            if len(info["SAMPLEGENO"]) != len(self.opt.add_genoinfo):
                self.raise_exception("SAMPLEGENO is unmatched.")
            else:
                print('\tOK. SAMPLEGENO='+','.join(info["SAMPLEGENO"]))

        except KeyError:
            self.raise_exception("No SAMPLEGENO.")

    def validate_variant_fields_in_only_vcf(self):
        self.load_mutanno_structure()

        for r1 in self.reader:
            # varkey = r1.CHROM +":"+ str(r1.POS) + "_" + r1.REF + ">" + ",".join(r1.ALT)
            self.varkey = r1.CHROM +":"+ str(r1.POS) 

            # print("checking..", self.varkey)
            for infokey in self.mutanno_sources.keys():
                if infokey in r1.INFO.keys():
                    flag = True
                    for attr in r1.INFO[infokey]:
                        
                        fields = attr.split('|')
                        # print(">>>>",fields)
                        if len(self.mutanno_sources[infokey]['fields']) != len(fields):
                            msg = infokey + " fields is not matched. (" + '|'.join(self.mutanno_sources[infokey]['fields']) + " != " + attr + ")" 
                            self.raise_exception(msg)

                        for v1 in fields:
                            if v1 != "":
                                flag = False
                    if flag:
                        msg = infokey + " has no fields ("+infokey+"=" + attr + ")" 
                        self.raise_exception(msg)                        
                            
                    self.check_vep_consequence_and_most_severe(r1.INFO, self.mutanno_sources['VEP']['fields'], r1)

            self.check_feature_ncbi(r1.INFO)
            self.check_genes(r1.INFO)
            if self.opt is not None and self.opt.add_genoinfo:
                self.check_samplegeno(r1.INFO)
            if self.opt is not None and self.opt.add_hg19:
                self.check_hg19(r1.INFO)

    def validate(self):
        self.metadata = read_metadata_and_check_duplicate(self.header)
        self.validate_variant_fields_in_only_vcf()
        # print(self.opt)
        # print(self.mutanno_sources)

def parse_tsi_annot_header(tsi_info_header):
    h = {}
    for f1 in tsi_info_header.split(';'):
        arr = f1.split('=')
        h[arr[0]] = arr[1].split('|')
    return h

def parse_info(infofield):
    d = {}
    for f1 in infofield.split(';'):
        if '=' in f1:
            arr = f1.split('=')
            d[arr[0]] = []
            for f2 in arr[1].split(','):
                d[arr[0]].append(f2.split('|'))
        else:
            if f1 != '':
                d[f1] = True
    return d

class AnnotTSIValidator(AnnotVCFValidator):

    def __init__(self, tsi_file):
        self.vcf_file = tsi_file
        self.header = ""
        self.colnames = []
        self.mutanno_sources = {}
        self.dslist = None
        self.opt = None
        self.varkey = ""
        self.annot_header = {}
        self.read_header()

    def read_header(self):
        for line in file_util.gzopen(self.vcf_file):
            line = file_util.decodeb(line)
            if line[0] != '#':
                break
            else:
                self.colnames = line.strip().split('\t')
                self.annot_header = parse_tsi_annot_header(self.colnames[-1])
                # print(self.annot_header)


    def validate_variant_fields_in_only_tsi(self):
        for line in file_util.gzopen(self.vcf_file):
            line = file_util.decodeb(line)
            if line[0] != '#':
                arr = line.split('\t')
                infofield = parse_info(arr[-1].strip())
                for k1 in infofield.keys():
                    # print(infofield[k1])
                    for d2 in infofield[k1]:
                        if len(self.annot_header[k1]) != len(d2):
                            msg = k1 + ": the field number not matched."
                            self.raise_exception(msg)

                # print('infofield:', infofield)
                self.check_vep_consequence_and_most_severe(infofield, self.annot_header['VEP'])
        

    def validate(self):
        self.validate_variant_fields_in_only_tsi()



