import tabix
import os
from ..util import file_util
from ..util import struct_util

def dv(dict1, key1, default):
    return struct_util.get_dict_value(dict1, key1, default)


class DataSource():
    def __init__(self, datastructure, sourcedir=""):
        self.datastructure = datastructure
        self.sourcefile = dv(self.datastructure, 'sourcefile', '')
        if sourcedir == "":
            self.sourcedir = dv(self.datastructure, 'sourcedir', '')
        else:
            self.sourcedir = sourcedir
        self.format = dv(self.datastructure, 'format', '')
        self.tchrom = ""
        self.header = ""
        self.colnames = []
        self.tb = None

    def load_datasourcefile(self, chrom):
        if self.tchrom != chrom:
            if "#CHROM#" in self.sourcefile or self.tb == None:
                self.sourcefile2 = os.path.join(self.sourcedir, self.sourcefile.replace("#CHROM#", chrom))
                self.tb = tabix.open(self.sourcefile2)
                self.tchrom = chrom
                self.load_header()
                self.set_format()

    def set_format(self):
        if self.format == "":
            if "|" in self.colnames[-1] and "CHROM" in self.colnames and "POS" in self.colnames:
                self.format = "tsi"

    def load_header(self):
        for line in file_util.gzopen(self.sourcefile2):
            line = file_util.decodeb(line)
            if line[0] == "#":
                self.header = line
                self.colnames = self.header[1:].split('\t')
            else:
                break

    def get_info(self, chrom, pos, ref, alt):
        info = ""
        rec = self.get_record(chrom, pos, ref, alt)
        if self.format == "tsi" and len(rec)>0:
            info = rec[-1]
        return info

    def get_record(self, chrom, pos, ref, alt):
        self.load_datasourcefile(chrom)
        rec = []
        try:
            recs = self.tb.query(chrom, pos - 1, pos + 1)
        except tabix.TabixError:
            recs = []
        for r1 in recs:
            if int(r1[1]) == pos and r1[3] == ref and r1[4] == alt:
                rec = r1
                break
        return rec


