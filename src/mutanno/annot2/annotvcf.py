import time
from ..util import file_util
from ..util import struct_util
# from ..util.struct_util import get_dict_value as dv
from ..util import vcf_util
from ..util import seq_util
from ..model.datasource import DataSourceList
from ..renderer import VCFRenderer
from ..renderer import JSONRenderer
from .. import _version

INFOIDX = vcf_util.VCF_COL.index('INFO')


class VCFAnnotator():
    vcfheader = None
    datastruct = {}
    opt = {}

    def __init__(self, opt):
        self.vcfheader = None

        self.opt = opt
        self.out = self.opt.out
        self.dslist = DataSourceList(self.opt.ds)
        if self.opt.sourcefile != "":
            self.dslist.sourcefile = self.opt.sourcefile
            self.dslist.sourcefile2 = self.opt.sourcefile
        if self.opt.single_source_mode:
            self.dslist.single_source_mode = self.opt.single_source_mode

        # self.add_hgvs = opt.add_hgvs
        # self.add_variant_class = opt['add_variant_class']
        # self.add_hg19 = opt['add_hg19']
        # self.add_genetable = opt['add_genetable']
        # self.split_multi_allelic_variant = opt['split_multi_allelic_variant']

        self.vcfreader = VCFReader(self.opt)
        self.vcfrenderer = VCFRenderer()
        self.jsonrenderer = JSONRenderer(self.dslist, self.opt.out)


    def get_datafile_header(self, datafile, sourcename, dataformat):
        global DATAHEADER
        header = []
        try:
            header = DATAHEADER[sourcename]
        except KeyError:
            for line in file_util.gzopen(datafile):
                line = line.decode('UTF-8')
                if dataformat == "tab":
                    if line[:len("#CHR")].lower() == "#chr":
                        arr = line[1:].strip().split('\t')
                        header = arr[5:]
                        break
                elif dataformat == "tsv":
                    if line[:len("#CHR")].lower() == "#chr":
                        arr = line[1:].strip().split('\t')
                        header = arr[4:]
                        break
                elif dataformat == "bed":
                    if line[:len("#CHR")].lower() == "#chr":
                        arr = line[1:].strip().split('\t')
                        header = arr
                        break
                elif dataformat == "info":
                    if line[:len("##INFO=")] == "##INFO=":
                        infoID = line[len("##INFO="):].split(',')[0].replace('<ID=', '')
                        target_fieldid = sourcename
                        if '.' in sourcename:
                            target_fieldid = sourcename.split('.')[1].strip()
                        if infoID == target_fieldid:
                            headercont = line.split(
                                'Description="')[-1].split(':')[-1].strip().replace('">', '').replace("'", "")
                            if ' | ' in headercont:
                                header = headercont.split(' | ')
                            else:
                                header = headercont.split('|')
                            break
                elif dataformat == "infovar":
                    if line[:len("##INFO=")] == "##INFO=":
                        infoID = line[len("##INFO="):].split(',')[0].replace('<ID=', '')
                        header.append(infoID)
                if line[0] != '#':
                    break
        DATAHEADER[sourcename] = header
        return header

    def get_annot_header(self):
        cont = ""


        for headertype in ["MUTANNO", "INFO"]:

            if self.opt.add_genoinfo != False:
                cont += vcf_util.get_info_header(headertype, "SAMPLEGENO", "v0.4", "05/07/2020" ,"Sample genotype information", ["NUMGT", "GT", "AD","SAMPLEID"], "samplegeno")
            
            mutanno_fields = []
            if self.opt.add_variant_class:
                mutanno_fields.append("variant_class")
            if self.opt.add_hgvs:
                mutanno_fields.append("hgvsg")
            if self.opt.split_multi_allelic_variant:
                # mutanno_fields.append("samplevariantkey")
                cont += vcf_util.get_info_header(headertype, "MULTIALLELE", _version.VERSION, _version.VERSION_DATE, "sample variant key for multiallelic variant", ["SAMPLEVARIANTKEY"])

            if len(mutanno_fields) > 0:
                cont += vcf_util.get_info_header(headertype, "MUTANNO", _version.VERSION, _version.VERSION_DATE, "mutanno", mutanno_fields)

            if self.opt.add_hg19:
                if self.opt.add_hgvs:
                    cont += vcf_util.get_info_header(headertype, "HG19", "v0.4", "05/07/2020","hg19 coordinates", ["chrom", "pos", "hgvsg"], "hg19")
                else:
                    cont += vcf_util.get_info_header(headertype, "HG19", "v0.4", "05/07/2020", "hg19 coordinates", ["chrom", "pos"], "hg19")
            if self.opt.add_genetable:
                cont += vcf_util.get_info_header(headertype, "GENES","v99","05/07/2020" ,"Gene table", ["ENSG","MOST_SEVERE_TRANSCRIPT","MOST_SEVERE_CONSEQUENCE"], "genes")

        if 'merged_one_field' in self.datastruct.keys() and self.datastruct['merged_one_field'] != '':
            fields = []
            for s1 in self.dslist.source_list:
                if s1.is_available:
                    if s1.fieldselection == "all":
                        if 'use_sourcename_in_fieldname' in self.datastruct.keys():
                            if not self.datastruct['use_sourcename_in_fieldname']:
                                fields.extend(DATAHEADER[s1['name']])
                        else:
                            for f1 in DATAHEADER[s1['name']]:
                                fields.append(s1['name'] + "_" + DATAHEADER[s1['name']][f1])

                    else:
                        for f1 in s1.field_list:
                            if 'use_sourcename_in_fieldname' in self.datastruct.keys() and \
                                    not self.datastruct['use_sourcename_in_fieldname']:
                                if 'name2' in f1.keys():
                                    fields.append(f1['name2'])
                                else:
                                    fields.append(f1['name'])
                            else:
                                if 'name2' in f1.keys():
                                    fields.append(s1['name'] + "_" + f1['name2'])
                                else:
                                    fields.append(s1['name'] + "_" + f1['name'])
            sourcename = self.datastruct['merged_one_field']
            sourcedesc = dv(self.datastruct, 'merged_one_field_desc', '')
            cont += vcf_util.get_info_header("INFO", sourcename, "","", sourcedesc, fields)
        else:
            for s1 in self.dslist.source_list:
                if s1.is_available:
                    fields = []
                    if s1.fieldselection == "all":
                        fields = DATAHEADER[s1.name]
                        fieldlist = []
                        for fieldname in fields:
                            fieldlist.append({'name': fieldname, 'desc': ''})
                        s1['fields'] = fieldlist
                    else:
                        for f1 in s1.field_list:
                            if f1.is_available:
                                fields.append(f1.name2)
                    cont += vcf_util.get_info_header("INFO",s1.name,"","", s1.desc, fields, s1.subembedded)
        return cont

    def get_version_info(self):
        d = {}
        d['MUTANNO'] = {}
        # d['MUTANNO']['MUTANNO'] = {'Version':_version.VERSION, 'Date':_version.VERSION_DATE}
        for s1 in self.dslist.source_list:
            if s1.is_available:
                d['MUTANNO'][s1.name] = {}
                if s1.version != "":
                    d['MUTANNO'][s1.name]['Version'] = s1.version
                if s1.version_date != "":
                    d['MUTANNO'][s1.name]['Date'] = s1.version_date
        return vcf_util.convert_to_metadata(d)

    def open_outpointer(self):
        self.out = file_util.strip_gzext(self.opt.out)
        file_util.check_dir(self.out)
        self.fp = open(self.out, 'w')

    def close_outpointer(self):
        self.fp.close()
        print("Saved " + self.out)


    def write_header(self):
        header = []
        headercont, colheader = self.vcfreader.get_header(keep_only_main_chrom=True)
        header.append(headercont.strip())
        header.append(self.get_version_info().strip())
        header.append(self.get_annot_header().strip())
        header.append(colheader.strip())
        self.fp.write('\n'.join(header) + '\n')

    def run(self):
        self.open_outpointer()
        stime0 = time.time()
        total_varno = 0
        self.write_header()

        ivar = 0
        stime1 = time.time()
        while(not self.vcfreader.eof):
            variant = self.vcfreader.get_variant()
            if self.vcfreader.eof or variant == None:
                break
            variant.set_annot(self.dslist)
            self.fp.write(self.vcfrenderer.render_vcfvariant(variant) + '\n')
            self.jsonrenderer.save_jsonvariant(variant)
            
            ivar += 1
            if ivar % 1000 == 0:
                etime1 = time.time()
                elapsed1 = etime1 - stime1
                log = 'processed... ' + str(ivar) + '/' + str(self.vcfreader.total_variant)
                log += " elapsed:" + str(round(elapsed1, 2)) + "s"
                stime1 = time.time()
                print(log)

        etime0 = time.time()
        elapsed0 = etime0 - stime0
        if ivar > 0:
            print("total:", elapsed0, ", time per variant:", elapsed0 / ivar,
                  ", time for 4M variants:", str(round(elapsed0 / ivar * 4000000 / 3600, 3)) + "hr")
        
        self.close_outpointer()
        
        if self.opt.out.endswith('.gz'):
            file_util.save_gzip(file_util.strip_gzext(self.opt.out))


class VCFReader():
    def __init__(self, opt={}):
        self.vcf = opt.vcf
        self.opt = opt 
        self.add_genoinfo = opt.add_genoinfo
        self.add_hgvs = opt.add_hgvs
        self.add_variant_class = opt.add_variant_class
        
        self.fp = file_util.gzopen(self.vcf)
        self.eof = False
        self.total_variant = 0
        self.i_variant = 0
        self.sampleid_list = []
        self.liftover_hg38_hg19 = None
        if self.opt.add_hg19:
            self.liftover_hg38_hg19 = seq_util.load_liftover(self.opt.chain)

    def is_clean_tag(self, line):
        flag = False
        for prevtag in ["##MUTANNO=<ID=", "##INFO=<ID="]:
            if line[:len(prevtag)] == prevtag:
                for tag in self.opt.clean_tag_list:
                    if tag in line:
                        flag = True
        return flag

    def get_header(self, keep_only_main_chrom=True):
        headercont = ""
        colheader = ""
        for line in file_util.gzopen(self.vcf):
            line = file_util.decodeb(line)
            if line[:2] == "##":
                if line[:len('##contig=')] == '##contig=' and keep_only_main_chrom:
                    contigchrom = line.split(',')[0].replace('##contig=<ID=chr','')
                    if len(contigchrom) < 3:
                        headercont += line
                else:
                    if not self.is_clean_tag(line):
                        headercont += line
            elif line[0] == "#":
                colheader = line
                self.colnames = line[1:].strip().split('\t')
                self.sampleid_list = self.colnames[9:]
            elif line.strip() != '':
                self.total_variant += 1
        return headercont, colheader

    def get_variant(self):
        variant = None
        while True:
            line = self.fp.readline()
            line = file_util.decodeb(line)

            if line.strip() != '':
                if line[0] != '#':
                    variant = VCFVariant(line, self.colnames, self.opt, self.liftover_hg38_hg19)
                    break
            else:
                self.eof = True
                break

        return variant

class VCFVariant:
    def __init__(self, vcfline="", colnames=[], opt={}, liftover_hg38_hg19=None, is_subvariant=False):
        self.opt = opt
        self.vcfline = vcfline
        self.record = vcfline.split('\t')
        self.record[-1] = self.record[-1].strip()
        self.record[INFOIDX] = vcf_util.strip_info(self.record[INFOIDX])
        self.vcfinfo = VCFVariantINFO(self.record[INFOIDX])
        self.chrom = self.record[0]
        self.nchrom = self.record[0].replace('chr','')
        self.pos = int(self.record[1])
        self.ref = self.record[3]
        self.alt = self.record[4]
        self.annotlines = []
        self.colnames = colnames
        self.sampleid_list = colnames[9:]
        self.is_multiallelic = False
        self.is_subvariant = is_subvariant
        self.variant_class = vcf_util.get_variant_class(self.ref, self.alt)
        self.split_variants = []
        self.samplevariantkey = self.chrom + ':' + str(self.pos) + '%20' + self.ref + '/' + self.alt
        self.annotdata = {}
        self.liftover_hg38_hg19 = liftover_hg38_hg19
        self.set_option()
        # print('processing...',self)

    def __str__(self):
        return self.chrom + ':' + str(self.pos)  + '_'  + self.ref + '>' + self.alt

    def set_option(self):
        try:
            if self.opt.add_genoinfo != False and not self.is_subvariant:
                self.add_genotypeinfo_in_infofield(self.opt.add_genoinfo)

            if len(self.opt.clean_tag_list) > 0:
                self.clean_tag_infofield(self.opt.clean_tag_list)

            if self.opt.add_variant_class or self.opt.add_hgvs or self.opt.split_multi_allelic_variant:
                self.add_mutanno_infofield(self.opt.add_variant_class, self.opt.add_hgvs, self.opt.split_multi_allelic_variant)

            if self.opt.add_hg19:
                self.add_hg19_in_infofield(self.liftover_hg38_hg19, self.opt.add_hgvs)
        except AttributeError:
            pass

    def add_genotypeinfo_in_infofield(self, target_sampleid):
        if len(target_sampleid) == 0:
            target_sampleidx = range(9, len(self.colnames))
        else:
            target_sampleidx = [self.sampleid_list.index(sid)+9 for sid in target_sampleid]

        samplegeno = []
        for k in target_sampleidx:
            gtinfo = self.record[k].split(':')
            converted_gtinfo = []
            converted_gtinfo.append(gtinfo[0].replace('|', '/'))
            converted_gtinfo.append(vcf_util.get_genotype(gtinfo[0], self.record[3], self.record[4]))
            converted_gtinfo.append(gtinfo[1].replace(',', '/'))
            converted_gtinfo.append(self.sampleid_list[k-9])
            samplegeno.append('|'.join(converted_gtinfo))
        self.record[INFOIDX] = vcf_util.add_info(self.record[INFOIDX], "SAMPLEGENO=" + ",".join(samplegeno))

    def add_mutanno_infofield(self, add_variant_class=False, add_hgvs=False, split_multi_allelic_variant=False):
        values = []
        if add_variant_class:
            values.append(vcf_util.encode_infovalue(self.variant_class))
        if add_hgvs:
            hgvsg = vcf_util.get_hgvsg(self.chrom, self.pos, self.ref, self.alt)
            values.append(vcf_util.encode_infovalue(hgvsg))
        if split_multi_allelic_variant:
            self.set_split_multiallelic_variant()
            # values.append(self.samplevariantkey)
        if len(values) > 0:
            self.record[INFOIDX] = vcf_util.add_info(self.record[INFOIDX], "MUTANNO=" + '|'.join(values))

        
    def add_hg19_in_infofield(self, liftover_hg38_hg19, add_hgvs=False):
        hg19 = ""
        for hg19coord in seq_util.convert_coordinate(liftover_hg38_hg19, self.chrom, self.pos):
            hg19_chrom = hg19coord[0]
            hg19_pos = hg19coord[1]
            if hg19 != "":
                hg19 += ","
            hg19 += hg19_chrom + '|' + str(hg19_pos)
            if add_hgvs:
                hg19 += '|' + vcf_util.get_hgvsg(hg19_chrom, hg19_pos, self.ref, self.alt)
        if hg19 != "":
            self.record[INFOIDX] = vcf_util.add_info(self.record[INFOIDX], "HG19=" + hg19)
    
    def clean_tag_infofield(self, clean_tag_list):
        rstinfo = ""
        for infofield in self.record[INFOIDX].split(';'):
            tag = infofield.split('=')[0]
            if tag not in clean_tag_list:
                if rstinfo != "":
                    rstinfo += ";"
                rstinfo += infofield
        self.record[INFOIDX] = rstinfo

    def set_split_multiallelic_variant(self):
        if "," in self.alt:
            for record in vcf_util.split_multiallelic_variants(self.record):
                record[INFOIDX] = vcf_util.add_info(record[INFOIDX], "MULTIALLELE=" + vcf_util.encode_infovalue(self.samplevariantkey))
                vcfline = '\t'.join(record)
                subvariant = VCFVariant(vcfline, self.colnames, self.opt, self.liftover_hg38_hg19, is_subvariant=True)

                self.split_variants.append(subvariant)
            self.is_multiallelic = True

    def set_annot(self, dslist):
        self.dslist = dslist
        if self.is_multiallelic:
            for variant in self.split_variants:
                variant.set_annot(dslist)
        else:
            self.annotmerger = dslist.get_variant_annot(self.nchrom, self.pos, ref=self.ref, alt=self.alt, variant=self)
        
    def get_vcf_info(self):
        self.vcfinfo.set_datasourcelist(self.dslist)
        return self.vcfinfo.get_info_dict()



class VCFVariantINFO:
    def __init__(self, recordstring=""):
        self.recordstring = recordstring
        self.data = None
        self.vcfheader = None
        self.dslist = None

    def set_datasourcelist(self, dslist):
        self.dslist = dslist

    def get_info_dict(self):
        if self.data is None:
            if self.dslist is None:
                self.parse_with_only_line()
            elif self.vcfheader is not None:
                self.parse_with_vcfheader()
            elif self.dslist is not None:
                self.parse_with_datasourcelist()
        return self.data

    def parse_with_vcfheader(self):
        pass

    def parse_with_datasourcelist(self):
        pass

    def parse_with_only_line(self):
        d = {}
        for field in self.recordstring.split(';'):
            if '=' in field:
                arr = field.strip().split('=')
                d[arr[0].strip()] = arr[1].split(',')
            else:
                d[field] = True
        self.data = d



