
from .util import vcf_util
from .util.vcf_util import INFOIDX
import json


class VCFRenderer:
    def __init__(self):
        pass

    def merge_twoinfofields(self, infofield1, infofield2):
        infofield1 = vcf_util.strip_info(infofield1)
        infofield2 = vcf_util.strip_info(infofield2)
        if infofield1 == "":
            info = infofield2
        if infofield2 == "":
            info = infofield1
        if infofield1 != "" and infofield2 != "":
            info = infofield1 + ';' + infofield2
        return info

    def convert_annotdata_to_infofield(self, annotdata, available_field_list):
        info = []
        for s1 in available_field_list.keys():
            infoattr = []
            for attr in annotdata[s1]:
                infofield = []
                for fieldname in available_field_list[s1]:
                    if isinstance(attr[fieldname], list):
                        for i in range(len(attr[fieldname])):
                            attr[fieldname][i] = vcf_util.encode_infovalue(attr[fieldname][i])
                        infofield.append("~".join(attr[fieldname]))
                    else:
                        attr[fieldname] = vcf_util.encode_infovalue(attr[fieldname])
                        infofield.append(attr[fieldname])
                infoattr.append('|'.join(infofield))
            if len(infoattr) > 0:
                merged_value = ','.join(infoattr)
                if merged_value.strip().replace("|", "") != "":
                    info.append(s1 + '=' + merged_value)
        return ';'.join(info)

    def render_vcfvariant(self, vcfvariant):
        if vcfvariant.is_multiallelic:
            lines = []
            for variant in vcfvariant.split_variants:
                record = variant.record
                info_record = self.render_vcfvariant_info(variant)
                record[INFOIDX] = self.merge_twoinfofields(record[INFOIDX], info_record)
                lines.append('\t'.join(record))
            return '\n'.join(lines)
        else:
            record = vcfvariant.record
            info_record = self.render_vcfvariant_info(vcfvariant)
            record[INFOIDX] = self.merge_twoinfofields(record[INFOIDX], info_record)
        return '\t'.join(record)

    def render_vcfvariant_info(self, vcfvariant):
        annotdata = vcfvariant.annotmerger.get_data()
        available_field_list = vcfvariant.annotmerger.available_field_list
        return self.convert_annotdata_to_infofield(annotdata, available_field_list)

    def render_vcfvariant_fast(self, vcfvariant, dslist):
        if vcfvariant.is_multiallelic:
            lines = []
            for variant in vcfvariant.split_variants:
                record = variant.record
                info_record = dslist.get_singlesource_annot(variant)
                record[INFOIDX] = self.merge_twoinfofields(record[INFOIDX], info_record)
                lines.append('\t'.join(record))
            return '\n'.join(lines)
        else:
            record = vcfvariant.record
            info_record = dslist.get_singlesource_annot(vcfvariant)
            record[INFOIDX] = self.merge_twoinfofields(record[INFOIDX], info_record)
        return '\t'.join(record)


class JSONRenderer:
    def __init__(self, dslist=None, outname=None):
        self.fp = None
        self.dslist = dslist
        self.i = 0
        if outname is not None:
            self.fp = open(outname + '.json', 'w')
            self.fp.write('[')

    def __del__(self):
        if self.fp is not None:
            self.fp.write(']')
            self.fp.close()

    def save_jsonvariant(self, vcfvariant):
        if self.i > 0:
            self.fp.write(',')
        self.fp.write(self.render_jsonvariant(vcfvariant))
        self.i += 1

    def render_jsonvariant(self, vcfvariant):
        if vcfvariant.is_multiallelic:
            records = []
            for variant in vcfvariant.split_variants:
                record = self.render_jsonvariant(variant)
                records.append(record)
            return ','.join(records)
        else:
            annotdata = vcfvariant.annotmerger.data
            annotdata['CHROM'] = vcfvariant.chrom
            annotdata['POS'] = vcfvariant.pos
            annotdata['REF'] = vcfvariant.ref
            annotdata['ALT'] = vcfvariant.alt
            rst = json.dumps(annotdata, indent=2)
        return rst


class TSIRenderer(VCFRenderer):
    def __init__(self):
        pass

    def render_vcfvariant(self, vcfvariant, keep_empty_variant=True, print_linereturn=True):
        annotdata = vcfvariant.annotmerger.get_data()
        available_field_list = vcfvariant.annotmerger.available_field_list
        info_record = self.convert_annotdata_to_infofield(annotdata, available_field_list)
        rst = ''
        if keep_empty_variant or info_record.strip() != '':
            cont = []
            cont.append(vcfvariant.chrom)
            cont.append(str(vcfvariant.pos))
            cont.append('')
            cont.append(vcfvariant.ref)
            cont.append(vcfvariant.alt)
            cont.append(info_record)
            rst = '\t'.join(cont)
            if print_linereturn:
                rst += '\n'
        return rst
