from ..util import file_util
from ..util import struct_util
from . import datastructure
import tabix


class MutannoDataParser:
    def __init__(self):
        pass

    def set_header(self, line):
        self.header = file_util.line2arr(line)
        self.header[0] = self.header[0].replace('#','')
        self.annotheader = {}
        for f1 in self.header[-1].split(';'):
            arrf1 = f1.split('=')
            self.annotheader[arrf1[0]] = arrf1[1].split('|')

    def parse(self, line):
        arr = file_util.line2arr(line)
        dat = {}
        for i in range(len(arr)-1):
            dat[self.header[i]] = arr[i]
        dat['POS'] = int(dat['POS'])
        for f1 in arr[-1].split(';'):
            arrf1 = f1.split('=')
            d = []

            for f11 in arrf1[1].split(','):
                arrf2 = f11.split('|')
                d2 = {}
                for i in range(len(self.annotheader[arrf1[0]])):
                    d2[self.annotheader[arrf1[0]][i]] = arrf2[i]
                d.append(d2)
            dat[arrf1[0]] = d
        return dat

class MutannoDataSequentialReader:
    def __init__(self, source_file):
        self.source_file = source_file
        self.fp = file_util.gzopen(source_file)
        self.mtparser = MutannoDataParser()
        self.load_header()

    def __iter__(self):
        return self

    def load_header(self):
        line = file_util.decodeb(self.fp.readline())
        if line[0] == '#':
            self.mtparser.set_header(line)

    def __next__(self):
        line = file_util.decodeb(self.fp.readline())
        if line.strip() == '':
            raise StopIteration
        elif line[0] == '#':
            self.mtparser.set_header(line)
        else:
            mtdata = self.mtparser.parse(line)
            return (mtdata)


class DataSource:
    def __init__(self):
        self.tabixpointer = ""
        self.chrom = ""
        self.name = ""
        self.name2 = ""
        self.ds = {}
        self.subheader = {}

    def set_datastructure(self, ds):
        self.ds = ds
        self.datapath = self.ds['datapath']
        self.name = self.ds['name']
        self.ref_column_index = struct_util.get_dict_value(self.ds, 'ref_column_index', 3)
        self.alt_column_index = struct_util.get_dict_value(self.ds, 'alt_column_index', 4)
        self.set_name2()

    def set_name2(self):
        self.name2 = struct_util.get_dict_value(self.ds, 'name2', "")
        if self.name2 == "":
            self.name2 = self.name

    def set_header(self, sourcefile):
        for line in file_util.gzopen(sourcefile):
            line = file_util.decodeb(line)
            if line[0] == "#":
                if line[:2] != "##":
                    self.header = file_util.line2arr(line[1:])
                    for j, h1 in enumerate(self.header):
                        if '|' in h1:
                            self.subheader[j] = h1.split('|')
                    break

    def tabix_load(self, chrom=""):
        sourcefile = self.datapath.replace('#CHROM#', chrom)
        if file_util.is_exist(sourcefile):
            print("loading..",sourcefile)
            self.tabixpointer = tabix.open(sourcefile)
            self.chrom = chrom
            self.set_header(sourcefile)

    def conv_data_to_dict(self, rec):
        d = struct_util.conv_dict_from_arr(rec, self.header)
        for j in self.subheader.keys():
            rec[j].split(',')
        
        rst.append(d)

    def get_sourcedata_variant(self, chrom, pos, ref, alt):
        rst = []
        if self.chrom != chrom:
            self.tabix_load(chrom)

        chrompos = chrom + ':' + str(pos) + '-' + str(pos)
        try:
            for rec in self.tabixpointer.querys(chrompos):
                if rec[self.ref_column_index] == ref and rec[self.alt_column_index] == alt:
                    d = self.conv_data_to_dict(rec)
                    rst.append(d)
        except AttributeError:
            pass
        except tabix.TabixError:
            pass

        return rst


class DataSourceList:
    def __init__(self):
        self.tabixpointer = {}
        self.variant_source_list = []

    def set_datastructure_json(self, dsfile):
        self.ds = datastructure.DataSourceListStructure(dsfile)

    def set_datastructure(self, dsobj):
        self.ds = dsobj

    def set_variant_source_list(self):
        if len(self.variant_source_list) == 0 and len(self.ds.variant_source_list) > 0:
            for s1 in self.ds.variant_source_list:
                datsrc = DataSource()
                datsrc.set_datastructure(s1)
                self.variant_source_list.append(datsrc)

    def get_sourcedata_variant(self, chrom, pos, ref, alt):
        self.set_variant_source_list()
        chrompos = chrom + ':' + str(pos) + '-' + str(pos)
        annot = {}
        for datsrc in self.variant_source_list:
            print(datsrc.name)
            recs = datsrc.get_sourcedata_variant(chrom, pos, ref, alt)

            # try:
            #     self.tabixpointer[s1['name']]
            # except KeyError:
            #     if file_util.is_exit(s1['datapath'].replace('#CHROM#', chrom)):
            #         self.tabixpointer[s1['name']] = tabix.open(s1['datapath'].replace('#CHROM#', chrom))

            # rec = self.tabixpointer[s1['name']].querys(chrompos)
            annot[datsrc.name] = recs

        return annot