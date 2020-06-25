from .._options import get_opt_object_from_dict
from ..util import file_util, vcf_util
import tabix


class TSIReaderSRC:
    def __init__(self, opt={}):
        if isinstance(opt, dict):
            self.opt = get_opt_object_from_dict(opt)
        else:
            self.opt = opt
        self.tsi = self.opt.tsi
        self.tb = tabix.open(self.tsi)
        self.tsiheader = []
        self.annotheader = {}
        self.headline = ""
        self.fp = file_util.gzopen(self.tsi)
        self.load_header()

    def load_header(self):
        for line in self.fp:
            line = file_util.decodeb(line)
            if line[0] == "#":
                self.headline = line
                self.tsiheader = line[1:].split('\t')
                self.tsiheader[-1] = self.tsiheader[-1].strip()
                self.annotheader = vcf_util.pars_info_header(self.tsiheader[-1])
                break

    def get_variant(self):
        line = file_util.decodeb(self.fp.readline())
        if line.strip() == "":
            v1 = None
        else:
            v1 = TSIVariant(line, self.tsiheader, self.annotheader)
        return v1

    def get_variants(self, genomic_range):
        varlist = []
        for rec in self.tb.querys(genomic_range):
            line = '\t'.join(rec) + '\n'
            if line.strip() == "":
                v1 = None
            else:
                v1 = TSIVariant(line, self.tsiheader, self.annotheader)
            varlist.append(v1)
        return varlist


class TSIVariant:
    def __init__(self, line, tsiheader, annotheader):
        record = line.split('\t')
        record[-1] = record[-1].strip()
        self.line = line

        self.chrom = record[0]
        self.pos = int(record[1])
        self.id = record[2].strip()
        self.ref = record[3].strip()
        self.alt = record[4].strip()
        self.qual = record[5]
        self.filter = record[6]
        self.annot = {}
        self.annotheader = annotheader

        self.is_range = False
        self.spos = self.pos
        self.epos = self.pos
        if 'END=' in record[-1]:
            for f1 in record[-1].split(';'):
                if f1[:len('END=')] == 'END=':
                    self.epos = int(f1.replace('END=', ''))
                    break
            self.is_range = True
        self.record = record
        self.load_annot(record[-1])

    def __str__(self):
        if self.is_range:
            s1 = self.chrom + ':' + str(self.spos) + '-' + str(self.epos)
            # if self.ref != "" and self.alt != "":
            #     s1 += '_' + self.ref + '/' + self.alt
        else:
            s1 = self.chrom + ':' + str(self.pos) + '_' + self.ref + '/' + self.alt
        return s1

    def load_annot(self, cont):
        d = {}
        for sourcefield in cont.split(";"):
            arr = sourcefield.split('=')
            sname = arr[0]
            if sname != 'END':
                try:
                    d[sname]
                except KeyError:
                    d[sname] = []
                for attr in arr[1].split(','):
                    d2 = {}
                    arr2 = attr.split('|')
                    # print(sname, len(arr2), len(self.annotheader[sname]), arr2, self.annotheader[sname])
                    try:
                        for idx, fieldname in enumerate(self.annotheader[sname]):
                            d2[fieldname] = arr2[idx].strip()
                    except IndexError:
                        print("Error: Field number is not matching." + str(self) + " " + sname +
                              " source. " + str(len(self.annotheader[sname])) + " " + str(len(arr2)))
                        print(self.annotheader[sname])
                        print(arr2)
                    d[sname].append(d2)
        self.annot = d

    def get_line(self):
        annotlist = []

        if self.ref == '':
            annotlist.append('END=' + str(self.epos))

        for k1 in self.annot.keys():
            attr2 = []
            for attr in self.annot[k1]:
                flist = []
                for f1 in self.annotheader[k1]:
                    flist.append(attr[f1])
                attr2.append('|'.join(flist))
            annotlist.append(k1 + '=' + ','.join(attr2))

        self.record[-1] = ';'.join(annotlist)
        line = '\t'.join(self.record) + '\n'
        return line
