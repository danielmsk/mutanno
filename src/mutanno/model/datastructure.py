from ..util import file_util
from ..util.struct_util import get_dict_value as dv
from .. import _version
from . import datasource
import os


# (ATTRIBUTE_NAME, DEFAULT_VALUE, IS_ESSENTIAL)
DSLIST_ATTR = [
    ("name", "MUTANNO", True),
    ("version", _version.VERSION, True),
    ("version_date", _version.VERSION_DATE, True),
    ("sourcefile", "", False),
    ("sourcefile_path", "", False),
    ("single_source_mode", False, False),
    ("level", "", True),
    ("source", [], True),
    ("hgvs", False, False),
    ("variant_class", False, False),
    ("hg19", False, False),
    ("fixpl", False, False),
    ("clean_tag_list", [], False),
]

DS_ATTR = [
    ("name", "", True),
    ("desc", "", True),
    ("version", "", True),
    ("version_date", "", True),
    ("sourcefile", "", True),
    ("function", None, True),
    ("param", None, True),
    ("subembedded", "", True),
    ("format", "tsi", True),
    ("chrompre", "", True),
    ("fieldselection", "", True),
    ("is_available", True, True),
    ("datatype", "variant", True),
    ("ref_column_index", 3, True),
    ("alt_column_index", 4, True),
    ("fields", [], True)
]

DSFIELD_ATTR = [
    ("name", "", True),
    ("name2", "", True),
    ("type", "string", True, ["string", "number", "integer", "boolean"]),
    ("desc", "", False),
    ("is_list", False, True),
    ("delimiter", "~", True),
    ("function", None, True),
    ("param", None, True),
    ("default", None, True),
    ("subembedded_key", False, True),
    ("keep_equal_str", [], True),
    ("is_available", True, True),
]


class DataSourceListStructure:
    def __init__(self, dsjsonfile=""):
        if isinstance(dsjsonfile, str) and file_util.is_exist(dsjsonfile):
            self.source_list = []
            self.available_source_list = []
            self.sources = {}
            self.dsjsonfile = dsjsonfile
            self.sourcefile2_list = []
            self.sourcelist_with_default = []
            self.load_dsfile()
            self.update_sourcefile2()
            self.fast_mapping_mode = False
            if len(self.sourcefile2_list) == 0:
                self.single_source_mode = True
                self.fast_mapping_mode = True

    def load_dsfile(self):
        ds = file_util.load_json(self.dsjsonfile)
        for attr in DSLIST_ATTR:
            if attr[0] == "source":
                for source_ds in ds["source"]:
                    dss = datasource.DataSource(source_ds)
                    dss.update_sourcefile2(self.sourcefile_path)
                    if dss.sourcefile2 != '':
                        self.sourcefile2_list.append(dss.sourcefile2)

                    # print(">>>>>>DSS:", dss, dss.is_available)
                    if dss.is_available:
                        self.source_list.append(dss)
                        self.available_source_list.append(dss)
                        self.sources[dss.name] = dss
                        if dss.has_default_value:
                            self.sourcelist_with_default.append(dss)
            else:
                self.__dict__[attr[0]] = dv(ds, attr[0], attr[1])

    def update_sourcefile2(self):
        if self.sourcefile_path != "":
            self.sourcefile2 = os.path.join(self.sourcefile_path, self.sourcefile)
        else:
            self.sourcefile2 = self.sourcefile


class DataSourceStructure:
    def __init__(self, dsjson=None):
        self.field_list = []
        self.available_field_list = []
        self.default_value_list = []
        self.has_default_value = False
        self.fields = {}
        self.fields2 = {}
        self.field_name2name1 = {}
        if isinstance(dsjson, dict):
            for attr in DS_ATTR:
                if attr[0] == "fields":
                    for field_ds in dsjson["fields"]:
                        dsfs = DataSourceFieldStructure(field_ds)
                        self.field_list.append(dsfs)
                        self.fields[dsfs.name] = dsfs
                        self.field_name2name1[dsfs.name2] = dsfs.name

                        if dsfs.is_available:
                            self.available_field_list.append(dsfs)
                            if dsfs.default is not None:
                                self.default_value_list.append(dsfs.default)
                                self.has_default_value = True
                            else:
                                self.default_value_list.append('')
                else:
                    self.__dict__[attr[0]] = dv(dsjson, attr[0], attr[1])

    def __str__(self):
        return self.name

    def update_sourcefile2(self, sourcefile_path):
        # print(">DataSourceStructure.update_sourcefile2()", self.name)
        if sourcefile_path != "":
            self.sourcefile2 = os.path.join(sourcefile_path, self.sourcefile)
        else:
            self.sourcefile2 = self.sourcefile
        # print(self.sourcefile2)


class DataSourceFieldStructure:
    def __init__(self, dsjson=None):
        if isinstance(dsjson, dict):
            for attr in DSFIELD_ATTR:
                self.__dict__[attr[0]] = dv(dsjson, attr[0], attr[1])
            self.set_name2()

    def __str__(self):
        return self.name

    def set_name2(self):
        if self.name2 == "":
            self.name2 = self.name

    def convert_type_from_string(self, strvalue):
        if self.type == "number":
            if strvalue == "":
                rst = self.default
            else:
                rst = float(strvalue)
        elif self.type == "integer":
            if strvalue == "":
                rst = self.default
            else:
                rst = int(strvalue)
        elif self.type == "boolean":
            # >> convert boolean value in external function. 
            # if strvalue == "1":
            #     rst = True
            # elif strvalue == "0":
            #     rst = False
            # else:
            #     rst = None
            rst = strvalue
        else:
            if strvalue == "":
                rst = self.default
            else:
                rst = strvalue
        return rst
