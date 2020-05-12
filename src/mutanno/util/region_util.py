

class Region:
    def __init__(self, region_str):
        self.region_str = region_str
        self.chrom = ""
        self.spos = 0
        self.epos = 0

    def pars_region_str():
        r = None
        if region != "":
            arr = region.split(':')
            arr2 = arr[1].split('-')
            r = {}
            self.chrom = arr[0].replace('chr', '')
            self.spos = int(arr2[0])
            self.epos = int(arr2[1])
        