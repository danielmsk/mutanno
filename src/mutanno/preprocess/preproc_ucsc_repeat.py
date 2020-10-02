import os
from ..util import file_util, vcf_util


class PreprocUCSCRepeat():
    def __init__(self, data):
        self.selected_fields = []

        if 'rawfiles' in data.keys():
            self.sqlfile = data['rawfiles'][0]
            self.infile = data['rawfiles'][1]
            self.dirpath = data['dirpath']
            self.outfile_title = data['outfile_title']
            self.out = os.path.join(self.dirpath, self.outfile_title)
        else:
            self.sqlfile = data['infile'][0]
            self.infile = data['infile'][1]
            self.out = data['out']
        self.set_outfile_extension()

    def set_outfile_extension(self):
        out2 = self.out.replace('.bed.gz', '.bed')
        if out2[-4:] != ".bed":
            out2 += ".bed"
        self.out = out2

    def get_header(self, fieldnames):
        cont = ['chrom','chromStart','chromEnd','name']
        for k in range(len(fieldnames)):
            if not fieldnames[k] in ['chrom','chromStart','chromEnd','name']:
                cont.append(fieldnames[k])
        return "#" + '\t'.join(cont)

    def load_fieldnames(self):
        fieldnames = []
        flag = False
        for line in open(self.sqlfile):		
            if line[:2] != '/*' and line[:2] != '--':
                line = line.strip()
                if line != "":
                    if "KEY" in line:
                        flag = False
                    if flag and line[0] == '`':
                        arr = line.replace('`','').split(' ')
                        # print (arr[0])
                        fieldnames.append(arr[0])
                    if "CREATE TABLE" in line:
                        flag = True

        return fieldnames

    def load_posname(self, fieldnames):
        m = {}
        if 'genoName' in fieldnames:
            m['chrom'] = "genoName"
            m['spos'] = "genoStart"
            m['epos'] = "genoEnd"
        else:
            m['chrom'] = "chrom"
            m['spos'] = "chromStart"
            m['epos'] = "chromEnd"
        return m

    def convert_to_bed(self):
        fieldnames = self.load_fieldnames()
        # print(fieldnames)
        posname = self.load_posname(fieldnames)
        
        f = open(self.out, 'w')
        f.write(self.get_header(fieldnames) + "\n")
        i = 0
        for line in file_util.gzopen(self.infile):
            line = file_util.decodeb(line)
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()

            m = {}
            for k in range(len(fieldnames)):
                m[fieldnames[k]] = arr[k].strip()
            cont = []
            cont.append(m[posname['chrom']].replace('chr',''))
            cont.append(m[posname['spos']])
            cont.append(m[posname['epos']])
            for k in range(len(fieldnames)):
                if not fieldnames[k] in ['chrom','chromStart','chromEnd']:
                    cont.append(m[fieldnames[k]])
            f.write('\t'.join(cont)+'\n')
        f.close()

        out = self.out
        print("Bzipping and Tabixing...", out)
        file_util.save_tabixgz(out)
        print("Saved...", out + ".gz")
        file_util.check_and_remove(out + '.gz.tbi', out, 3)


    def run(self):
        self.convert_to_bed()
        
