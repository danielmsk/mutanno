from . import seq_util


class Region:
    def __init__(self, region_str=""):
        self.region_str = region_str
        self.chrom = ""
        self.spos = 0
        self.epos = 0
        self.block_spos = 0
        self.block_epos = 0
        self.has_chr = False
        self.pars_region_str()

    def __str__(self):
        if self.region_str == "all":
            return "all_region"
        else:
            return self.chrom + ':' + str(self.spos) + '-' + str(self.epos)

    def pars_region_str(self):
        if self.region_str == "all":
            self.spos = 1
            self.epos = seq_util.CHROM_LEN['hg38']['1']
            self.chrom = '1'
        elif self.region_str != "":
            arr = self.region_str.split(':')
            arr2 = arr[1].split('-')
            if 'chr' in arr[0]:
                self.has_chr = True
            self.chrom = arr[0]
            self.nchrom = arr[0].replace('chr', '')
            self.spos = int(arr2[0])
            self.epos = int(arr2[1])

    def get_split_region(self, bin_size=10000):
        regions = []
        for k in range(self.spos, self.epos, bin_size):
            r1 = Region()
            r1.chrom = self.chrom
            r1.spos = k
            r1.epos = r1.spos + bin_size - 1
            if r1.epos > self.epos:
                r1.epos = self.epos
            regions.append(r1)
        return regions

    def get_next_block_region(self, block_size):
        if self.block_spos == 0:
            self.block_spos = self.spos
        else:
            self.block_spos = self.block_epos + 1

        if self.block_spos > self.epos:
            if self.region_str == "all":
                if self.chrom == "M":
                    next_block = None
                else:
                    chromlist = seq_util.CHROM_LEN['hg38'].keys()
                    self.chrom = chromlist[chromlist.index(self.chrom) + 1]
                    self.spos = 1
                    self.epos = seq_util.CHROM_LEN['hg38'][self.chrom]
            else:
                next_block = None
        else:
            if (self.block_spos + block_size) <= self.epos:
                self.block_epos = self.block_spos + block_size - 1
            else:
                self.block_epos = self.epos
            next_block = Region(self.chrom + ':' + str(self.block_spos) + '-' + str(self.block_epos))
        return next_block
