
from .util import region_util
from .util import vcf_util
from .util import file_util
from .annotvcf import VCFReader, VCFVariant
import tabix
import json


class ANNOTVCFViewer(VCFReader):
    def __init__(self, opt):
        self.opt = opt
        self.vcf = opt.vcf
        self.region = region_util.Region(opt.region)
        print(opt.region)

        # inherited value
        self.infoheader = {}
        self.clean_tag_list = []
        self.total_variant = 0
        self.tb = tabix.open(self.vcf)

    def get_variant(self, region):
        variants = []
        for recs in self.tb.querys(region.region_str):
            variant = VCFVariant(recs, self.colnames, self.opt, False, infoheader=self.infoheader)
            variants.append(variant)
            # print(variant.vcfinfo.get_info_dict())
        return variants

    def run(self):
        headercont, colheader = self.get_header()
        self.colnames = colheader.split('\t')

        variants = self.get_variant(self.region)
        for v1 in variants:
            pass
            # print(v1.vcfinfo.infoheader['GNOMAD'])
            print(json.dumps(v1.vcfinfo.get_info_dict(), sort_keys=False, indent=2, separators=(',', ': ')))

        #

        # print (self.vcf)
        # print(self.region)
