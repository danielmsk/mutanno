
from ..util import file_util
from ..util import struct_util
import os



class DataSourceListStructure:
    def __init__(self, dsjsonfile=""):
        self.ds = {}
        self.variant_source_list = []
        if dsjsonfile != "":
            self.ds = file_util.load_json(dsjsonfile)
        self.load_sourcefile()

    def load_sourcefile(self):
        datafile_path = struct_util.get_dict_value(self.ds, 'datafile_path', "")

        for s1 in self.ds['source']:
            datafile = struct_util.get_dict_value(s1, 'datafile', "")
            s1['datapath'] = os.path.join(datafile_path, datafile)
            self.variant_source_list.append(s1)

