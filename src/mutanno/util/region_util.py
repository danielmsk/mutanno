

class Region:
    def __init__(self, region_str=""):
        self.region_str = region_str
        self.chrom = ""
        self.spos = 0
        self.epos = 0
        self.pars_region_str()

    def __str__(self):
        return self.chrom + ':' + str(self.spos) + '-' + str(self.epos)

    def pars_region_str(self):
        r = None
        if self.region_str != "":
            arr = self.region_str.split(':')
            arr2 = arr[1].split('-')
            r = {}
            self.chrom = arr[0].replace('chr', '')
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
