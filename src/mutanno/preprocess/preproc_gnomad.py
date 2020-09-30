import tabix
import os
from ..util import file_util, seq_util, vcf_util

class PreprocGnomAD():
    def __init__(self, data):
        self.header = ""
        self.fieldlist = []
        self.dirpath = data['dirpath']
        self.outfile_title = data['outfile_title']
        self.rawfiles = data['rawfiles']
        self.refversion = data['refversion']
        self.source_name = data['source_name']
        self.out = self.get_out_filename()
        self.rawfileset = {}
        self.set_rawfiles_by_chrom()

    def get_out_filename(self):
        out = os.path.join(self.dirpath, self.outfile_title + ".chr#CHROM#" + '.mti')
        return out

    def set_header(self):
        chromlist = list(self.rawfileset.keys())
        tb = tabix.open(self.rawfileset[chromlist[0]])
        i = 0
        fieldset = {}
        for rec in tb.querys("chr"+chromlist[0] + ":1-10000000"):
            i += 1
            for f1 in rec[7].split(';'):
                if '=' in f1 and f1[:len('vep=')] != 'vep=':
                    arr = f1.split('=')
                    fieldset[arr[0]] = 1
            if i > 1000:
                break
        self.fieldlist = list(fieldset.keys())
        arrheader = ["CHROM", "POS", "ID", "REF", "ALT"]
        arrheader.append(self.source_name + "=" + "|".join(self.fieldlist))
        self.header = "#" + '\t'.join(arrheader)

    def get_header(self):
        if self.header == "":
            self.set_header()
        return self.header

    def convert_mti_record(self, rec):
        mti_rec = []
        mti_rec.append(rec[0].replace('chr', ''))
        mti_rec.append(rec[1])
        mti_rec.append(rec[2])
        mti_rec.append(rec[3])
        mti_rec.append(rec[4])

        d = {}
        for f1 in rec[7].split(';'):
            if "=" in f1:
                arr = f1.split('=')
                d[arr[0]] = vcf_util.encode_value(arr[1])

        datlist = []
        for fn in self.fieldlist:
            try:
                datlist.append(d[fn])
            except KeyError:
                datlist.append("")
        mti_rec.append('|'.join(datlist))
        return mti_rec

    def convert_to_mti_chrom(self, chrom):
        out = self.out.replace('#CHROM#', chrom)
        f = open(out, "w")
        f.write(self.get_header() + "\n")
        tb = tabix.open(self.rawfileset[chrom])
        chromlen = seq_util.CHROM_LEN[self.refversion][chrom]
        i = 0
        for rec in tb.querys("chr"+chrom + ":1-" + str(chromlen)):
            mti_rec = self.convert_mti_record(rec)
            f.write('\t'.join(mti_rec) + '\n')
            if i % 100000 == 0:
                print("Processing...", rec[0] + ':' + rec[1])
            i += 1
        f.close()
        
        print("Bzipping and Tabixing...", out)
        file_util.save_tabixgz(out)
        print("Saved...", out + ".gz")
        file_util.check_and_remove(out + '.gz.tbi', out, 3)

    def convert_to_mti(self):
        for chrom in seq_util.MAIN_CHROM_LIST:
            if chrom in self.rawfileset.keys():
                self.convert_to_mti_chrom(chrom)
        
    def set_rawfiles_by_chrom(self):
        for rawfile in self.rawfiles:
            if rawfile.endswith('.vcf.bgz'):
                for namefield in rawfile.split('/')[-1].split('.'):
                    if "chr" in namefield:
                        chrom = namefield.replace('chr','')
                        self.rawfileset[chrom] = rawfile

    def run(self):
        self.convert_to_mti()
        
