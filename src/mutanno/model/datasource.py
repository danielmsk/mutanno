from os.path import isfile
from ..util import file_util
from ..util import struct_util
from . import datastructure
from .. import external_functions
from .datastructure import DataSourceStructure, DataSourceListStructure
import tabix


def run_external_function(function, params, paramvalues):
    # print("paramvalues", paramvalues['mutanno_value_sections'])

    exec_str = "external_functions." + function
    exec_str += '('
    paramstr_list = []
    for param in params:
        # print('param:', param)
        if param in paramvalues.keys():
            paramstr_list.append('paramvalues["' + param + '"]')
            # print('paramstr_list:', paramstr_list)
        else:
            paramstr_list.append(param)
    exec_str += ','.join(paramstr_list)
    exec_str += ')'
    # print('@@@@',exec_str)
    # print(paramvalues)
    rstvalue = eval(exec_str)
    return rstvalue


class VariantAnnotMerger:
    def __init__(self):
        self.data = {}
        self.annot = {}
        self.annot_list = []
        self.available_field_list = {}
        self.flag_set_data = False

    def add_variant_annot(self, variant_annot):
        self.annot[variant_annot.source.name] = variant_annot
        self.annot_list.append(variant_annot)

    def get_data(self):
        if not self.flag_set_data:
            self.set_data()
        return self.data

    def set_data(self):
        for annot in self.annot_list:
            annotdata = annot.get_data(self.data)
            sname = annot.source.name
            try:
                self.data[sname]
            except KeyError:
                self.data[sname] = []
            self.data[sname].extend(annotdata)

            try:
                self.available_field_list[sname].extend(annot.available_field_list)
            except KeyError:
                self.available_field_list[sname] = annot.available_field_list
        self.flag_set_data = True


class VariantAnnot:
    def __init__(self, records, header_record, datasource=None, return_all=False, variant=None) :
        self.data = []
        self.rawdata = []
        self.records = records
        self.header_record = header_record
        self.header = []
        self.source = datasource # DataSource
        self.return_all = return_all
        self.available_field_list = []
        self.flag_set_data = False
        self.othersource_annotdata = {}
        self.variant = variant
        self.vcf_info = None
        if self.variant is not None:
            self.vcf_info = self.variant.get_vcf_info()
        self.set_rawdata()

    def get_data(self, othersource_annotdata):
        self.othersource_annotdata = othersource_annotdata
        if not self.flag_set_data:
            self.set_data()
        return self.data

    def set_rawdata(self):
        # print("### set_rawdata():" + self.source.name)
        if len(self.rawdata) == 0:
            if self.source.format == "tsi":
                self.set_rawdata_from_tsi()
            if self.source.format == "bed":
                if self.return_all or self.source.is_available:
                    self.set_rawdata_from_bed()
            if self.source.format == "tsv":
                if self.return_all or self.source.is_available:
                    self.set_rawdata_from_tsv()
            if self.source.format == "function":
                if self.return_all or self.source.is_available:
                    self.set_rawdata_from_function()

    def set_data(self):
        self.set_data_from_rawdata()
        self.flag_set_data = True

    def set_data_from_rawdata(self):
        # print("### set_data_from_rawdata(): " + self.source.name)
        if self.source.is_available:
            for f1 in self.source.available_field_list:
                self.available_field_list.append(f1.name)
            # print(self.rawdata)
            # print("\t===(1)self.rawdata:", self.rawdata)
            self.rawdata = self.apply_external_function(self.rawdata)
            # print("\t===(2)self.rawdata:", self.rawdata)
            for attr in self.rawdata:
                d = {}
                for f1 in self.source.available_field_list:
                    # if self.source.fields[f1.name].function is None:
                    try:
                        d[f1.name] = attr[f1.name]
                        # print(">>>", f1.name, attr[f1.name])
                    except KeyError:
                        if f1.default != None:
                            d[f1.name] = f1.default
                        else:  
                            d[f1.name] = ""
                self.data.append(d)
            if len(self.data) == 0:
                d = {}
                for f1 in self.source.available_field_list:
                    if f1.default != None:
                        d[f1.name] = f1.default
                if len(d.keys()) > 0:
                    self.data = [d]
            # print("\tself.data:", self.data)


    def apply_external_function(self, sections):
        # print(">>>sections",sections)
        for section_idx, section in enumerate(sections):
            # print(">>>>section", section)

            if self.source.function is not None:
                function = self.source.function
                params = self.source.param.split(',')
                paramvalues = {
                    "mutanno_value_sections":sections,
                    "mutanno_value_section_idx":section_idx,
                    "mutanno_value_columnheader":self.header,
                    "vcf_info_value":self.vcf_info,
                    "mutanno_value_data":self.othersource_annotdata
                }
                for k1 in section.keys():
                    paramvalues[k1] = section[k1]
                sections = run_external_function(function, params, paramvalues)

            for field in self.source.available_field_list:
                try:
                    if self.source.fields[field.name].function is not None:

                        field = self.source.fields[field.name]
                        function = field.function
                        # print("########", function)
                        params = field.param.split(',')
                        paramvalues = {
                            "mutanno_value_sections":sections,
                            "mutanno_value_section_idx":section_idx,
                            "mutanno_value_columnheader":self.header,
                            "vcf_info_value":self.vcf_info,
                            "mutanno_value_data":self.othersource_annotdata
                        }
                        for k1 in section.keys():
                            paramvalues[k1] = section[k1]
                        section[field.name] = run_external_function(function, params, paramvalues)
                except KeyError:
                    pass
        return sections
    
    def set_rawdata_from_tsi(self):
        """
        The raw data is unselected values from data source.
        """
        if '=' in self.header_record[-1]:
            for s1 in self.header_record[-1].split(';'):
                arr = s1.split('=')
                source_name = arr[0]
                if source_name == self.source.name:
                    values = arr[1].split('|')
                    self.header = arr[1].split('|')
        else:
            self.header = self.header_record[-1].split('|')

        self.header = self.source.conv_fieldname1(self.header)

        for rec in self.records:
            for s1 in rec[-1].split(';'):
                if '=' in s1:
                    arr = s1.split('=')
                    source_name = arr[0]
                    r1s = arr[1].split(',')
                    source_is_available = True
                else:
                    source_name = self.source.name
                    source_is_available = self.source.is_available
                    r1s = s1.split(',')

                if source_name == self.source.name:
                    
                    if self.return_all or source_is_available:
                        for r1 in r1s:
                            values = r1.split('|')

                            d = self.set_init_annot_section_key(rec)
                            for i, fieldname in enumerate(self.header):
                                d[fieldname] = self.convert_field_type(fieldname, values[i])
                                
                            if len(d.keys()) > 0:
                                self.rawdata.append(d)


    def set_init_annot_section_key(self, rec):
        d = {}
        d['CHROM'] = rec[0]
        d['POS'] = int(rec[1])
        d['REF'] = rec[self.source.ref_column_index]
        d['ALT'] = rec[self.source.alt_column_index]
        d['VARKEY'] = d['CHROM'] + ':' + str(d['POS']) + '_' + d['REF'] + '>' + d['ALT']
        return d

    def check_field_availability(self, fieldname):
        flag = False
        try:
            self.source.fields[fieldname]
            if self.return_all or self.source.fields[fieldname].is_available:
                flag = True
        except KeyError:
            flag = self.return_all
        return flag

    def convert_field_type(self, fieldname, value):
        try:
            fieldstructure = self.source.fields[fieldname]
            if fieldstructure.is_list:
                if value == "":
                    rst = None
                else:
                    rst = []
                    for v1 in value.split(fieldstructure.delimiter):
                        rst.append(fieldstructure.convert_type_from_string(v1))
            else:
                rst = fieldstructure.convert_type_from_string(value)
        except KeyError:
            if value == "":
                rst = None
            else:
                rst = value
        return rst

    def set_rawdata_from_bed(self):
        # print("???1", self.source.name, self.records)
        self.header = []
        for rec in self.records:
            d = self.set_init_annot_section_key(rec)
            for i, fieldname in enumerate(self.header_record):
                # if self.check_field_availability(fieldname):
                d[fieldname] = self.convert_field_type(fieldname, rec[i])
                self.header.append(fieldname)
            self.rawdata.append(d)

    def set_rawdata_from_tsv(self):
        self.header = []
        for rec in self.records:
            d = self.set_init_annot_section_key(rec)
            for i, fieldname in enumerate(self.header_record):
                # if self.check_field_availability(fieldname):
                d[self.header_record[i]] = self.convert_field_type(fieldname, rec[i])
                self.header.append(fieldname)
            self.rawdata.append(d)

    def set_rawdata_from_function(self):
        # print(">VariantAnnot.set_rawdata_from_function() : " + self.source.name)
        d = {}
        for field in self.source.field_list:
            d[field.name] = "tmpOOOOOOOOOO"
            self.header.append(field.name)
        self.rawdata.append(d)

    
class DataSource(DataSourceStructure):
    tchrom = ""
    source_name = ""
    header = ""
    tabixpointer = None
    
    def set_header(self, sourcefile):
        # print(">DataSource.set_header():", self.name)
        if self.header == "":
            for line in file_util.gzopen(sourcefile):
                line = file_util.decodeb(line)
                if line[0] == "#":
                    if line[:2] != "##":
                        self.header = file_util.line2arr(line[1:])
                        break

    def conv_fieldname1(self, tmp_name_list):
        if isinstance(tmp_name_list, str):
            tmp_name_list = [tmp_name_list]
        rstlist = []
        for tname in tmp_name_list:
            try:
                rstlist.append(self.field_name2name1[tname])
            except KeyError:
                rstlist.append(tname)
        return rstlist

    def tabix_load(self, chrom=""):
        # print(">DataSource.tabix_load():", self.name)
        sourcefile = self.sourcefile2.replace('#CHROM#', chrom)
        # print(">>", sourcefile, self.name)
        if file_util.is_exist(sourcefile) and isfile(sourcefile):
            # print("loading..", sourcefile)
            self.tabixpointer = tabix.open(sourcefile)
            self.tchrom = chrom
            self.set_header(sourcefile)
        else:
            if self.sourcefile != "":
                raise Exception("Source file is not exist. : " + sourcefile )

    def get_variant_annot(self, chrom, pos, ref="", alt="", epos=-9, return_raw=True, records = [], header_record = [], variant=None):
        # print(self.sourcefile2)
        if self.tchrom != chrom and (self.tabixpointer == None  or "#CHROM#" in self.sourcefile2):
            self.tabix_load(chrom)

        if len(header_record) > 0:
            recs = records
            header = header_record
        else:
            if epos > pos:
                chrompos = chrom + ':' + str(pos) + '-' + str(epos)
            else:
                chrompos = chrom + ':' + str(pos) + '-' + str(pos)
            recs = []
            try:
                for rec in self.tabixpointer.querys(chrompos):
                    if (epos < 0 and rec[self.ref_column_index] == ref and rec[self.alt_column_index] == alt
                        or (epos > pos)
                    ):
                        recs.append(rec)
            except AttributeError:
                pass
            except tabix.TabixError:
                pass
            header = self.header

        return VariantAnnot(recs, header, datasource=self, variant=variant)

    def get_variantkeys(self, chrom, spos, epos):
        if self.tchrom != chrom and (self.tabixpointer == None  or "#CHROM#" in self.sourcefile2):
            self.tabix_load(chrom)
        chrompos = chrom + ':' + str(spos) + '-' + str(epos)
        vkeys = {}
        if self.tabixpointer is not None:
            for rec in self.tabixpointer.querys(chrompos):
                vkey = rec[0] + "_" + rec[1] + "_" + "_" + rec[self.ref_column_index] + "_" + rec[self.alt_column_index]
                pos = int(rec[1])
                try:
                    vkeys[pos]
                except KeyError:
                    vkeys[pos] = []
                if vkey not in vkeys[pos]:
                    vkeys[pos].append(vkey)
            
        return vkeys



class DataSourceList(DataSourceListStructure):

    def get_variant_annot(self, chrom, pos, ref="", alt="", epos=-9, variant=None):
        chrompos = chrom + ':' + str(pos) + '-' + str(pos)
        annotmerger = VariantAnnotMerger()

        if self.sourcefile != "": 
            singlesource_ds = {
                "name": "MUTANNO_SINGLESOURCE",
                "sourcefile": self.sourcefile,
                "format": "tsi",
                "is_available": True,
                "fields": []
            }
            self.single_source = DataSource(singlesource_ds)
            self.single_source.sourcefile2 = self.sourcefile2
            single_source_annot = self.single_source.get_variant_annot(chrom, pos, ref, alt, return_raw=True)
            # annotmerger.add_variant_annot(single_source_annot)

        for s1 in self.source_list:
            if self.single_source_mode:
                annot = s1.get_variant_annot(chrom, pos, ref=ref, alt=alt, epos=epos, records=single_source_annot.records, header_record=single_source_annot.header_record, variant=variant )
            else:
                annot = s1.get_variant_annot(chrom, pos, ref=ref, alt=alt, epos=epos, variant=variant)
            annotmerger.add_variant_annot(annot)

        return annotmerger

    def get_variantkeys(self, chrom, spos, epos):
        merged_vkeys = {}
        for s1 in self.source_list:
            if s1.format in ["tsi", "tsv"]:
                vkeys = s1.get_variantkeys(chrom, spos, epos)
                for pos in vkeys.keys():
                    try:
                        merged_vkeys[pos]
                    except KeyError:
                        merged_vkeys[pos] = []

                    for vkey in vkeys[pos]:
                        if vkey not in merged_vkeys[pos]:
                            merged_vkeys[pos].append(vkey)
        return merged_vkeys


# class MutannoDataParser:
#     def __init__(self):
#         pass

#     def set_header(self, line):
#         self.header = file_util.line2arr(line)
#         self.header[0] = self.header[0].replace('#','')
#         self.annotheader = {}
#         for f1 in self.header[-1].split(';'):
#             arrf1 = f1.split('=')
#             self.annotheader[arrf1[0]] = arrf1[1].split('|')

#     def parse(self, line):
#         arr = file_util.line2arr(line)
#         dat = {}
#         for i in range(len(arr)-1):
#             dat[self.header[i]] = arr[i]
#         dat['POS'] = int(dat['POS'])
#         for f1 in arr[-1].split(';'):
#             arrf1 = f1.split('=')
#             d = []

#             for f11 in arrf1[1].split(','):
#                 arrf2 = f11.split('|')
#                 d2 = {}
#                 for i in range(len(self.annotheader[arrf1[0]])):
#                     d2[self.annotheader[arrf1[0]][i]] = arrf2[i]
#                 d.append(d2)
#             dat[arrf1[0]] = d
#         return dat

# class MutannoDataSequentialReader:
#     def __init__(self, source_file):
#         self.source_file = source_file
#         self.fp = file_util.gzopen(source_file)
#         self.mtparser = MutannoDataParser()
#         self.load_header()

#     def __iter__(self):
#         return self

#     def load_header(self):
#         line = file_util.decodeb(self.fp.readline())
#         if line[0] == '#':
#             self.mtparser.set_header(line)

#     def __next__(self):
#         line = file_util.decodeb(self.fp.readline())
#         if line.strip() == '':
#             raise StopIteration
#         elif line[0] == '#':
#             self.mtparser.set_header(line)
#         else:
#             mtdata = self.mtparser.parse(line)
#             return (mtdata)

