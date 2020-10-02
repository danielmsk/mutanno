import os
from ..util import file_util, proc_util


class PreprocCADD():
    def __init__(self, data):
        self.selected_fields = []

        if 'rawfiles' in data.keys():
            self.infiles = data['rawfiles']
            self.dirpath = data['dirpath']
            self.outfile_title = data['outfile_title']
            self.out = os.path.join(self.dirpath, self.outfile_title)
        else:
            self.infile = data['infile'][0]
            self.out = data['out']

    def move_tsv(self):
        for infile in self.infiles:
            cmd = "mv " + self.infile + " " + self.dirpath
            proc_util.run_cmd(cmd)

    def run(self):
        # self.convert_to_mti()
        self.move_tsv()
        
