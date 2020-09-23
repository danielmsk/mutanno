import tabix
import os
from ..util import file_util, seq_util

class PreprocGnomAD():
    def __init__(self, data):
        self.dirpath = data['dirpath']
        self.outfile_title = data['outfile_title']
        self.out = self.get_out_filename()
        print(data)
        seq_util.MAIN_CHROM_LIST

    def get_out_filename(self):
        out = os.path.join(self.dirpath, self.outfile_title + '.mti')
        return out


    def run(self):
        print(self.data)
        open(self.data[0], "r")
