
from ..util import file_util, vcf_util

vcf_fields = ["ALLELE", "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL"]
mti_fields = ["ALLELE", "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL", "MAXDS"]


class PreprocSpliceAI():
    def __init__(self, data):
        self.selected_fields = []
        self.infile = data['infile']
        self.out = data['out']
        self.set_outfile_extension()
    
    def set_outfile_extension(self):
        out2 = self.out.replace('.mti.gz', '.mti')
        if out2[-4:] != ".mti":
            out2 += ".mti"
        self.out = out2

    def get_header(self):
        cont = ["CHROM","POS","ID","REF","ALT"]
        cont.append("SpliceAI=" + "|".join(mti_fields))
        return "#" + '\t'.join(cont)

    def cal_maxds(self, data):
        maxds = 0
        for sc_name in ["DS_AG", "DS_AL", "DS_DG", "DS_DL"]:
            ds = float(data[vcf_fields.index(sc_name)])
            if maxds < ds:
                maxds = ds
        return maxds

    def convert_to_mti(self):
        f = open(self.out, 'w')
        f.write(self.get_header() + "\n")
        i = 0
        for line in file_util.gzopen(self.infile):
            line = file_util.decodeb(line)
            if line[0] != "#":
                arr = line.split('\t')
                data = arr[7].strip().replace('SpliceAI=', '').split('|')
                data.append(str(self.cal_maxds(data)))
                arr2 = [arr[0], arr[1], arr[2], arr[3], arr[4]]
                arr2.append('|'.join(data))
                f.write('\t'.join(arr2) + '\n')

                if i % 100000 == 0:
                    print("Processing...", arr[0] + ':' + arr[1])
                i += 1

        f.close()

        out = self.out
        print("Bzipping and Tabixing...", out)
        file_util.save_tabixgz(out)
        print("Saved...", out + ".gz")
        file_util.check_and_remove(out + '.gz.tbi', out, 3)


    def run(self):
        self.convert_to_mti()
        
