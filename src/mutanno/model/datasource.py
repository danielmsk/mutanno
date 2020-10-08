#########################################################################
#
# Libraries
#
#########################################################################
from os.path import isfile
from ..util import file_util, vcf_util
from .datastructure import DataSourceStructure, DataSourceListStructure
from .. import external_functions
import tabix

#########################################################################
#
# Functions
#
#########################################################################
def run_external_function(function, params, paramvalues):
    ''' '''
    exec_str = "external_functions." + function
    exec_str += '('
    paramstr_list = []
    for param in params:
        if param in paramvalues.keys():
            paramstr_list.append('paramvalues["' + param + '"]')
        else:
            paramstr_list.append(param)
        #end if
    #end for
    exec_str += ','.join(paramstr_list)
    exec_str += ')'
    rstvalue = eval(exec_str)
    return rstvalue
#end def run_external_function

#########################################################################
#
# Objects
#
#########################################################################
#########################################################################
#
# VariantAnnotMerger
#       store VariantAnnot objects
#       methods to iterate through stored objects
#
#########################################################################
class VariantAnnotMerger:

    def __init__(self):
        self.data = {}
        self.annot = {}
        self.annot_list = []
        self.available_field_list = {}
        self.flag_set_data = False
    #end def __init__

    def add_variant_annot(self, variant_annot):
        ''' add VariantAnnot object '''
        self.annot[variant_annot.source.name] = variant_annot
        self.annot_list.append(variant_annot)
    #end def add_variant_annot

    def get_data(self):
        ''' '''
        if not self.flag_set_data:
            self.set_data()
        #end if
        return self.data
    #end def get_data

    def set_data(self):
        ''' '''
        for annot in self.annot_list:
            annotdata = annot.get_data(self.data) # here is where extract a specific source,
                                                  # annotdata contains only annotation for
                                                  # source associated to annot
            ## annot is a VariantAnnot object
            ##  -> VariantAnnot.get_data()
            ##    -> VariantAnnot.set_data()
            ##       use source information from DataSource object associated to VariantAnnot
            ##       to retrieve only specific fields for source from rawdata that contains
            ##       all annotations for all sources
            sname = annot.source.name
            try:
                self.data[sname]
            except KeyError:
                self.data[sname] = []
            #end try
            self.data[sname].extend(annotdata)

            try:
                self.available_field_list[sname].extend(annot.available_field_list)
            except KeyError:
                self.available_field_list[sname] = annot.available_field_list
            #end try
        #end for
        self.flag_set_data = True
    #end def set_data

#end class VariantAnnotMerger

#########################################################################
#
# VariantAnnot
#       store annotations for a source
#       store source information
#
#########################################################################
class VariantAnnot:

    def __init__(self, records, header_record, datasource=None, return_all=False, variant=None, tmp_clinvar_idmap={}):
        self.data = []
        self.rawdata = []
        self.records = records
        self.header_record = header_record
        self.header = []
        self.source = datasource  # DataSource
        self.return_all = return_all
        self.available_field_list = []
        self.flag_set_data = False
        self.othersource_annotdata = {}
        self.variant = variant
        self.vcf_info = None
        if self.variant is not None:
            self.vcf_info = self.variant.get_vcf_info()
        #end if
        self.set_rawdata()
        self.tmp_clinvar_idmap = tmp_clinvar_idmap
    #end def __init__

    def get_data(self, othersource_annotdata):
        ''' '''
        self.othersource_annotdata = othersource_annotdata
        if not self.flag_set_data:
            self.set_data()
        #end if
        return self.data
    #end def get_data

    def set_rawdata(self):
        ''' '''
        if len(self.rawdata) == 0:
            if self.source.format == "tsi":
                self.set_rawdata_from_mti()
            #end if
            if self.source.format == "bed":
                if self.return_all or self.source.is_available:
                    self.set_rawdata_from_bed()
                #end if
            #end if
            if self.source.format == "tsv":
                if self.return_all or self.source.is_available:
                    self.set_rawdata_from_tsv()
                #end if
            #end if
            if self.source.format == "function":
                if self.return_all or self.source.is_available:
                    self.set_rawdata_from_function()
                #end if
            #end if
        #end if
    #end def set_rawdata

    def set_data(self):
        ''' '''
        self.set_data_from_rawdata()
        self.flag_set_data = True
    #end def set_data

    def set_data_from_rawdata(self):
        ''' '''
        # here is where extract the specific source annotations
        # source is DataSource object and contains
        # source name and source fields
        if self.source.is_available:
            for f1 in self.source.available_field_list:
                self.available_field_list.append(f1.name)
                # source fields from DataSource object are saved in this
                # VariantAnnot object available_field_list attribute
            #end for
            self.rawdata = self.apply_external_function(self.rawdata)
            for attr in self.rawdata: # raw data contains everything from data source file
                d = {}
                for f1 in self.source.available_field_list: # picking only fields for source
                    try:
                        d[f1.name] = attr[f1.name]
                    except KeyError:
                        if f1.default is not None:
                            d[f1.name] = f1.default
                        else:
                            d[f1.name] = ""
                        #end if
                    #end try
                #end for
                self.data.append(d)
            #end for
            if len(self.data) == 0:
                d = {}
                for f1 in self.source.available_field_list:
                    if f1.default is not None:
                        d[f1.name] = f1.default
                    #end if
                #end for
                if len(d.keys()) > 0:
                    self.data = [d]
                #end if
            #end if
        #end if
    #end def set_data_from_rawdata

    def apply_external_function(self, sections):
        for section_idx, section in enumerate(sections):
            if self.source.function is not None:
                function = self.source.function
                params = file_util.trim_split(self.source.param, ',')
                paramvalues = {
                    "mutanno_value_sections": sections,
                    "mutanno_value_section_idx": section_idx,
                    "mutanno_value_columnheader": self.header,
                    "vcf_info_value": self.vcf_info,
                    "tmp_clinvar_idmap": self.tmp_clinvar_idmap,
                    "mutanno_value_data": self.othersource_annotdata
                }
                for k1 in section.keys():
                    paramvalues[k1] = section[k1]
                sections = run_external_function(function, params, paramvalues)

            tmp_section = {}
            for field in self.source.available_field_list:
                try:
                    if self.source.fields[field.name].function is not None:
                        field = self.source.fields[field.name]
                        function = field.function
                        params = field.param.split(',')
                        paramvalues = {
                            "mutanno_value_sections": sections,
                            "mutanno_value_section_idx": section_idx,
                            "mutanno_value_columnheader": self.header,
                            "vcf_info_value": self.vcf_info,
                            "tmp_clinvar_idmap": self.tmp_clinvar_idmap,
                            "mutanno_value_data": self.othersource_annotdata
                        }
                        for k1 in section.keys():
                            paramvalues[k1] = section[k1]
                        tmp_section[field.name] = run_external_function(function, params, paramvalues)
                except KeyError:
                    pass

            for fieldname in tmp_section.keys():
                section[fieldname] = tmp_section[fieldname]

        return sections

    def set_rawdata_from_mti(self):
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
                    valuecont = arr[1].strip()
                    r1s = valuecont.split(',')
                    source_is_available = True
                else:
                    source_name = self.source.name
                    source_is_available = self.source.is_available
                    valuecont = s1.strip()
                    r1s = s1.split(',')

                if valuecont != '':
                    if source_name == self.source.name:
                        if self.return_all or source_is_available:
                            for r1 in r1s:
                                values = r1.split('|')
                                d = self.set_init_annot_section_key(rec)

                                # temporary code
                                if source_name == "CLINVAR_SUBMISSION":
                                    self.header = ['ClinVarAccession', 'Interpretation', 'DateLastEvaluated', 'ReviewStatus',
                                                   'AssertionMethod', 'AssertionMethodCitation', 'Method', 'Condition', 'ConditionXRef',
                                                   'AlleleOrigin', 'Submitter', 'SubmitterID', 'Citation', 'Comment']

                                for i, fieldname in enumerate(self.header):
                                    try:
                                        d[fieldname] = self.convert_field_type(fieldname, values[i])
                                    except IndexError:
                                        pass

                                if len(d.keys()) > 0:
                                    self.rawdata.append(d)

    def set_init_annot_section_key(self, rec):
        ''' '''
        d = {}
        d['CHROM'] = rec[0]
        d['POS'] = int(rec[1])
        # TODO:
        # d['REF'] = rec[self.source.ref_column_index]
        # d['ALT'] = rec[self.source.alt_column_index]
        # d['VARKEY'] = d['CHROM'] + ':' + str(d['POS']) + '_' + d['REF'] + '>' + d['ALT']
        return d
    #end def set_init_annot_section_key

    def check_field_availability(self, fieldname):
        ''' '''
        flag = False
        try:
            self.source.fields[fieldname]
            if self.return_all or self.source.fields[fieldname].is_available:
                flag = True
            #end if
        except KeyError:
            flag = self.return_all
        #end try
        return flag
    #end def check_field_availability

    def convert_field_type(self, fieldname, value):
        ''' '''
        try:
            fieldstructure = self.source.fields[fieldname]
            if fieldstructure.is_list:
                if value == "":
                    rst = None
                else:
                    rst = []
                    for v1 in value.split(fieldstructure.delimiter):
                        rst.append(fieldstructure.convert_type_from_string(v1))
                    #end for
                #end if
            else:
                rst = fieldstructure.convert_type_from_string(value)
            #end if
        except KeyError:
            if value == "":
                rst = None
            else:
                rst = value
            #end if
        return rst
    #end def convert_field_type

    def set_rawdata_from_bed(self):
        ''' '''
        self.header = []
        for rec in self.records:
            d = self.set_init_annot_section_key(rec)
            for i, fieldname in enumerate(self.header_record):
                # if self.check_field_availability(fieldname):
                d[fieldname] = self.convert_field_type(fieldname, rec[i])
                self.header.append(fieldname)
            #end for
            self.rawdata.append(d)
        #end for
    #end def set_rawdata_from_bed

    def set_rawdata_from_tsv(self):
        ''' '''
        self.header = []
        for rec in self.records:
            d = self.set_init_annot_section_key(rec)
            for i, fieldname in enumerate(self.header_record):
                # if self.check_field_availability(fieldname):
                d[self.header_record[i]] = self.convert_field_type(fieldname, rec[i])
                self.header.append(fieldname)
            #end for
            self.rawdata.append(d)
        #end for
    #end def set_rawdata_from_tsv

    def set_rawdata_from_function(self):
        ''' '''
        d = {}
        for field in self.source.field_list:
            d[field.name] = "tmpOOOOOOOOOO"
            self.header.append(field.name)
        #end for
        self.rawdata.append(d)
    #end def set_rawdata_from_function

#end class VariantAnnot

#########################################################################
#
# DataSource
#       object to store information for source
#       inherit from DataSourceStructure
#
#########################################################################
class DataSource(DataSourceStructure):

    tchrom = {}
    name = ""
    source_name = ""
    sourcefile = []
    sourcefile2 = []
    header = ""
    tabixpointer = {}
    ref_column_index = 3
    alt_column_index = 4

    def set_header(self, sourcefile):
        ''' '''
        if self.header == "":
            for line in file_util.gzopen(sourcefile):
                line = file_util.decodeb(line)
                if line[0] == "#":
                    if line[:2] != "##":
                        self.header = file_util.line2arr(line[1:])
                        self.infoheader = vcf_util.pars_info_header(self.header[-1])
                        break
                    #end if
                #end if
            #end for
        #end if
    #end def set_header

    def conv_fieldname1(self, tmp_name_list):
        ''' '''
        if isinstance(tmp_name_list, str):
            tmp_name_list = [tmp_name_list]
        #end if
        rstlist = []
        for tname in tmp_name_list:
            try:
                rstlist.append(self.field_name2name1[tname])
            except KeyError:
                rstlist.append(tname)
            #end try
        #end for
        return rstlist
    #end def conv_fieldname1

    def tabix_check_and_load(self, chrom=""):
        ''' open buffers to datasource files available as sourcefile '''
        for i, sf2 in enumerate(self.sourcefile2):
            if (i not in self.tchrom.keys() or self.tchrom[i] != chrom) and (i not in self.tabixpointer.keys() or self.tabixpointer[i] is None or "#CHROM#" in sf2):
                sourcefile = sf2.replace('#CHROM#', chrom)
                if file_util.is_exist(sourcefile) and isfile(sourcefile):
                    self.tabixpointer[i] = tabix.open(sourcefile)
                    self.tchrom[i] = chrom
                    self.set_header(sourcefile)
                else:
                    if self.sourcefile != "":
                        print("Source file is not exist. : " + sourcefile)
                    #end if
                #end if
            #end if
        #end for
    #end def tabix_check_and_load

    def get_variant_annot(self, chrom, pos, ref="", alt="", epos=-9, return_raw=True, records=[],
                          header_record=[], variant=None, tmp_clinvar_idmap={}):
        ''' access datasource files with tabix and get all annotations for variant '''
        ## !! this function is called for every source in available_source_list
        ## by DataSourceList.get_variant_annot(),
        ## every iteration save all annotations for all sources in records,
        ## created VariantAnnot objects have the same redundant records
        ## attribute that always contains annotations for all the sources !!
        self.tabix_check_and_load(chrom)
        recs = []
        for i in self.tabixpointer.keys():
            if len(header_record) > 0:
                recs = records
                header = header_record
            else:
                if epos > pos:
                    chrompos = chrom + ':' + str(pos) + '-' + str(epos)
                else:
                    chrompos = chrom + ':' + str(pos) + '-' + str(pos)
                #end if

                try:
                    for rec in self.tabixpointer[i].querys(chrompos):
                        ## !! this if is misterious, check why bed format is special !!
                        if ((epos < 0 and rec[self.ref_column_index] == ref and rec[self.alt_column_index] == alt)
                            or (epos > pos) or (self.format == "bed")
                                or (rec[self.ref_column_index] == "" and rec[self.alt_column_index] == "")):
                            # recs.append(rec)
                            if not rec in recs: recs.append(rec) # this should quick fix double annotations bug
                                                                 # when multiple datasource are used
                            #end if
                        #end if
                    #end for
                except AttributeError:
                    pass
                except tabix.TabixError:
                    pass
                #end try
                header = self.header
            #end if
        #end for
        return VariantAnnot(recs, header, datasource=self, variant=variant, tmp_clinvar_idmap=tmp_clinvar_idmap)
    #end def get_variant_annot

    def get_variantkeys(self, chrom, spos, epos):
        ''' '''
        self.tabix_check_and_load(chrom)
        chrompos = chrom + ':' + str(spos) + '-' + str(epos)
        vkeys = {}
        for i in self.tabixpointer.keys():
            try:
                for rec in self.tabixpointer[i].querys(chrompos):
                    p1 = int(rec[1])
                    if p1 >= spos and p1 <= epos:
                        vkey = rec[0] + "_" + rec[1] + "_" + "_" + \
                            rec[self.ref_column_index] + "_" + rec[self.alt_column_index]
                        pos = int(rec[1])
                        try:
                            vkeys[pos]
                        except KeyError:
                            vkeys[pos] = []
                        #end try
                        if vkey not in vkeys[pos]:
                            vkeys[pos].append(vkey)
                        #end if
                    #end if
                #end for
            except tabix.TabixError:
                pass
            #end try
        #end for
        return vkeys
    #end def get_variantkeys

#end class DataSource

#########################################################################
#
# DataSourceList
#       this is the main object the program works on
#       store all the information in the ds json file
#
#       calls several objects that update information at different levels:
#           -> DataSourceListStructure
#               -> DataSource
#                   -> DataSourceStructure
#                       -> DataSourceFieldStructure
#
#########################################################################
class DataSourceList(DataSourceListStructure):

    tabixpointer = {}
    tabixchrom = ""
    sourcefile = []
    ref_column_index = 3
    alt_column_index = 4
    merged_source = None
    use_raw_source = False
    merged_sourcefile_mode = False
    tmp_clinvar_idmap = {}  ## temporary

    def load_merged_source(self):
        ''' '''
        singlesource_ds = {
            "name": "MUTANNO_SINGLESOURCE",
            "sourcefile": self.sourcefile,
            "format": "tsi",
            "is_available": True,
            "fields": []
        }
        self.merged_source = DataSource(singlesource_ds)
        self.merged_source.sourcefile2 = self.sourcefile2
    #end def load_merged_source

    def get_variant_annot(self, chrom, pos, ref="", alt="", epos=-9, variant=None):
        ''' this create VariantAnnotMerger
        and populate it with VariantAnnot objects '''
        annotmerger = VariantAnnotMerger()
        if len(self.sourcefile) > 0:
            ## !! why we have to create and use a fake DataSource when we have the real one !!
            if self.merged_source is None:
                self.load_merged_source()
            #end if
            ## !! this is DataSource method that is being called
            ## same name is confusing !!
            merged_source_annot = self.merged_source.get_variant_annot(
                chrom, pos, ref, alt, return_raw=True, tmp_clinvar_idmap=self.tmp_clinvar_idmap)
            self.merged_sourcefile_mode = True
            # annotmerger.add_variant_annot(merged_source_annot)
        #end if
        for s1 in self.available_source_list:
            if self.merged_sourcefile_mode:
                annot = s1.get_variant_annot(chrom, pos, ref=ref, alt=alt, epos=epos,
                                             records=merged_source_annot.records,
                                             header_record=merged_source_annot.header_record, variant=variant, tmp_clinvar_idmap=self.tmp_clinvar_idmap)
                ## !! this call is just creating a new object
                ## with the same records from merged_source_annot
                ## it just change the source information in the new object !!
            else:
                annot = s1.get_variant_annot(chrom, pos, ref=ref, alt=alt, epos=epos,
                                             variant=variant, tmp_clinvar_idmap=self.tmp_clinvar_idmap)
            #end if
            annotmerger.add_variant_annot(annot)
        #end for
        return annotmerger
    #end def get_variant_annot

    def get_singlesource_annot(self, vcfvariant):
        ''' '''
        infos = []
        chrompos = vcfvariant.nchrom + ':' + str(vcfvariant.spos) + '-' + str(vcfvariant.epos)
        ref = vcfvariant.ref
        alt = vcfvariant.alt

        flag_annot = False
        for i in range(len(self.sourcefile)):
            try:
                tp = self.tabixpointer[i]
            except KeyError:
                tp = None
            #end try
            if tp is None or ("#CHROM#" in self.sourcefile[i] and self.tabixchrom != vcfvariant.nchrom):
                self.tabixpointer[i] = tabix.open(self.sourcefile[i].replace("#CHROM#", vcfvariant.nchrom))
                self.tabixchrom = vcfvariant.nchrom
            #end if
            try:
                for rec in self.tabixpointer[i].querys(chrompos):
                    if (((rec[self.ref_column_index] == ref and rec[self.alt_column_index] == alt)
                        or (rec[self.ref_column_index] == "" and rec[self.alt_column_index] == ""))
                            and vcfvariant.pos == int(rec[1])):
                        infos.append(rec[-1])
                        flag_annot = True
                    #end if
                #end for
            except AttributeError:
                pass
            except tabix.TabixError:
                pass
            #end try
            if flag_annot:
                break
            #end if
        #end for

        rst = ";".join(infos)
        for ds in self.sourcelist_with_default:
            if ds.name + '=' not in rst:
                infos.append(ds.name + '=' + '|'.join(ds.default_value_list))
                rst = ";".join(infos)
            #end if
        #end for
        return rst
    #end def get_singlesource_annot

    def get_variantkeys(self, chrom, spos, epos):
        ''' '''
        merged_vkeys = {}
        for s1 in self.source_list:
            if s1.format in ["tsi", "tsv"]:
                vkeys = s1.get_variantkeys(chrom, spos, epos)
                for pos in vkeys.keys():
                    try:
                        merged_vkeys[pos]
                    except KeyError:
                        merged_vkeys[pos] = []
                    #end try

                    for vkey in vkeys[pos]:
                        if vkey not in merged_vkeys[pos]:
                            merged_vkeys[pos].append(vkey)
                        #end if
                    #end for
                #end for
            #end if
        #end for
        return merged_vkeys
    #end def get_variantkeys

    def check_fast_mapping_mode(self):
        ''' '''
        fast_mapping_mode = False
        if self.use_raw_source:
            if self.merged_source is None:
                self.load_merged_source()
                self.merged_source.set_header(self.merged_source.sourcefile2[0].replace('#CHROM#', "1"))
            #end if

            fast_mapping_mode = True
            for s1 in self.available_source_list:
                if len(s1.available_field_list) == len(self.merged_source.infoheader[s1.name]):
                    for idx, f1 in enumerate(s1.available_field_list):
                        if f1.name != self.merged_source.infoheader[s1.name][idx]:
                            fast_mapping_mode = False
                            break
                        #end if
                    #end for
                    if not fast_mapping_mode:
                        break
                    #end if
                else:
                    fast_mapping_mode = False
                    break
                #end if
            #end for
        #end if
        self.fast_mapping_mode = fast_mapping_mode
        return fast_mapping_mode
    #end def check_fast_mapping_mode

#end class DataSourceList
