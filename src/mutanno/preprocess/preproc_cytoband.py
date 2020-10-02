import os
from ..util import file_util, vcf_util


class PreprocCytoband():
    def __init__(self, data):
        self.selected_fields = []

        if 'rawfiles' in data.keys():
            self.infile = data['rawfiles'][0]
            self.dirpath = data['dirpath']
            self.outfile_title = data['outfile_title']
            self.out = os.path.join(self.dirpath, self.outfile_title)
        else:
            self.infile = data['infile']
            self.out = data['out']
        self.set_outfile_extension()

    def set_outfile_extension(self):
        out2 = self.out.replace('.bed.gz', '.bed')
        if out2[-4:] != ".bed":
            out2 += ".bed"
        self.out = out2

    def get_header(self):
        cont = ["CHROM","SPOS","EPOS","cytoband"]
        return "#" + '\t'.join(cont)

    def convert_to_bed(self):
        f = open(self.out, 'w')
        f.write(self.get_header() + "\n")
        i = 0
        for line in file_util.gzopen(self.infile):
            line = file_util.decodeb(line)
            if line[0] != "#":
                arr = line.split('\t')
                arr2 = [arr[0].replace('chr', ''), arr[1], arr[2], arr[3]]
                f.write('\t'.join(arr2) + '\n')
        f.close()

        out = self.out
        print("Bzipping and Tabixing...", out)
        file_util.save_tabixgz(out)
        print("Saved...", out + ".gz")
        file_util.check_and_remove(out + '.gz.tbi', out, 3)


    def run(self):
        self.convert_to_bed()
        
