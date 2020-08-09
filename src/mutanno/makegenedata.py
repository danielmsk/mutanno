import os
import json
from .util import file_util
from .util import struct_util
from .util import seq_util
from . import external_functions

RESERVED_COL = ["chrom", "spos", "epos", "ensgid"]


def is_available(field):
    flag = True
    if struct_util.get_dict_value(field, "is_available", True) is False:
        flag = False
    if struct_util.get_dict_value(field, "do_import", True) is False:
        flag = False
    return flag


def is_reserved_column(f1):
    flag = False
    if f1['name'] in RESERVED_COL or ('name2' in f1.keys() and f1['name2'] in RESERVED_COL):
        flag = True
    return flag


def strip_value(v1):
    v1 = v1.replace('"', '')
    return v1


def encode_value(v1, delimiter=""):
    if v1 == '.' or v1 == '-':
        v1 = ''
    if delimiter != "":
        if delimiter != '|':
            v1 = v1.replace("|", "%7C")
            v1 = v1.replace(delimiter, '|')
    # v1 = urllib.parse.quote(v1)
    return v1


class GeneSourceReader():
    def __init__(self, datastruct, datafile_path):
        self.datastruct = datastruct
        self.source_name = self.datastruct['name']
        self.fname = os.path.join(datafile_path, self.datastruct['datafile'])
        self.fileformat = self.datastruct['format']
        self.target_colidx = []
        self.header = []
        self.subembedded = struct_util.get_dict_value(datastruct, 'subembedded', '')
        self.delimiter = {}
        self.list_type_fields = {}
        self.defaultvalue = {}
        self.field_function = {}
        self.field_function_param = {}
        self.field_names = []
        self.key_colidx = -999
        self.reserved_colidx = {}
        self.reserved_data = {}
        self.field_type = {}
        self.data = {}
        self.ensgid_list = []
        self.is_range_data = False
        self.target_colnames = []
        self.target_colnames2 = []
        self.is_ensemblegene = struct_util.get_dict_value(datastruct, 'is_ensemblegene', False)
        self.read_header()
        print("loading..", self.source_name)
        # print(self.header)
        self.set_target_column(self.datastruct['fields'])

        if self.key_colidx == -999:
            self.load_range_data()
            self.is_range_data = True
        else:
            self.load_data()

    def split_into_fields(self, line):
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        return arr

    def split_into_fields2(self, line):
        splited_fields = []
        for field1 in line.split('\t'):
            for f1 in field1.strip().split('|'):
                splited_fields.append(f1)
        return splited_fields

    def read_header(self):
        for line in file_util.gzopen(self.fname):
            if self.fname.endswith('.gz'):
                line = line.decode('UTF-8')
            if line[0] == '#' and line[1] != '#':
                self.header = self.split_into_fields(line[1:])
                break

    def set_target_column(self, fields_structure):
        name2_map = {}
        for f1 in fields_structure:

            if is_available(f1):
                self.field_names.append(f1['name'])

                # print("self.header", self.header, f1['name'], self.key_colidx)
                colidx = self.header.index(f1['name'])
                if f1['name'] in self.header:
                    if is_reserved_column(f1):
                        for resv in RESERVED_COL:
                            if f1['name'] == resv or ('name2' in f1.keys() and f1['name2'] == resv):
                                self.reserved_colidx[resv] = self.header.index(f1['name'])
                    else:
                        if 'is_list' in f1.keys() and f1['is_list']:
                            self.list_type_fields[self.header.index(f1['name'])] = f1['is_list']
                        if 'delimiter' in f1.keys():
                            self.delimiter[self.header.index(f1['name'])] = f1['delimiter']
                        if 'default' in f1.keys():
                            self.defaultvalue[self.header.index(f1['name'])] = f1['default']

                        self.target_colidx.append(colidx)
                else:
                    if not is_reserved_column(f1):
                        self.target_colidx.append(-999)

                if 'type' in f1.keys():
                    self.field_type[colidx] = f1['type']

                if 'function' in f1.keys() and f1['function'] != '':
                    self.field_function[colidx] = f1['function']
                    self.field_function_param[colidx] = ''
                    if 'param' in f1.keys() and f1['param'] != '':
                        self.field_function_param[colidx] = f1['param']

                name2_map[colidx] = struct_util.get_dict_value(f1, 'name2', '')

        if 'ensgid' in self.reserved_colidx.keys():
            self.key_colidx = self.reserved_colidx['ensgid']

        self.target_colnames = []
        self.target_colnames2 = []
        for tidx in self.target_colidx:
            self.target_colnames.append(self.header[tidx])
            if name2_map[tidx] == "":
                self.target_colnames2.append(self.header[tidx])
            else:
                self.target_colnames2.append(name2_map[tidx])

    def rm_version_in_ensgid(self, ensg):
        return ensg.strip().split('.')[0]

    def get_value_using_eval_field_function(self, colidx, cidx_no, cidxvalue):
        external_functions.void()
        exec_str = "external_functions." + self.field_function[colidx]
        exec_str += '('
        paramstr_list = []
        for param in self.field_function_param[colidx].split(','):
            if param == "mutanno_value_variantkey":
                pstr = 'variantkey'
            elif cidx_no == self.field_names.index(param):
                pstr = 'cidxvalue'
            else:
                pstr = 'arr_selected_fields['
                pstr += str(self.field_names.index(param))
                pstr += ']'
            paramstr_list.append(pstr)
        exec_str += ','.join(paramstr_list)
        exec_str += ')'
        cidxvalue = eval(exec_str)
        return cidxvalue

    def typecast(self, value, value_type):
        if (value == '' or value == 'NA') and value_type != 'string':
            rst_value = None
        else:
            if value_type == 'integer':
                try:
                    rst_value = int(value)
                except ValueError:
                    print('ValueError:', value)
                rst_value = int(value)
            elif value_type == 'number':
                try:
                    rst_value = float(value)
                except ValueError:
                    print('ValueError:', value)
                rst_value = float(value)
            elif value_type == 'boolean':
                rst_value = bool(value)
            else:
                rst_value = value
        return rst_value

    def load_data(self):
        for line in file_util.gzopen(self.fname):
            if self.fname.endswith('.gz'):
                line = line.decode('UTF-8')
            if line[0] != '#':
                arr = self.split_into_fields(line)

                ensgid = self.rm_version_in_ensgid(arr[self.key_colidx])
                self.ensgid_list.append(ensgid)
                try:
                    self.data[ensgid]
                except KeyError:
                    self.data[ensgid] = []

                cidx_no = 0
                target_data = []
                for colidx in self.target_colidx:
                    cidx_no += 1
                    if colidx == -999:
                        cidxvalue = ""
                    else:
                        cidxvalue = strip_value(arr[colidx])

                    if colidx in self.field_function.keys():
                        cidxvalue = self.get_value_using_eval_field_function(colidx, cidx_no, cidxvalue)

                    delimiter = ''
                    if colidx in self.delimiter.keys():
                        delimiter = self.delimiter[colidx]
                    if colidx in self.list_type_fields.keys():
                        cidxvaluelist = []
                        for v1 in cidxvalue.split(delimiter):
                            cidxvaluelist.append(self.typecast(v1, self.field_type[colidx]))
                        target_data.append(cidxvaluelist)
                    else:
                        target_data.append(self.typecast(cidxvalue, self.field_type[colidx]))

                self.data[ensgid].append(target_data)

                for resv in self.reserved_colidx.keys():
                    colidx = self.reserved_colidx[resv]
                    try:
                        self.reserved_data[resv]
                    except KeyError:
                        self.reserved_data[resv] = {}

                    try:
                        self.reserved_data[resv][ensgid]
                    except KeyError:
                        if resv == 'ensgid':
                            self.reserved_data[resv][ensgid] = self.rm_version_in_ensgid(arr[colidx])
                        else:
                            self.reserved_data[resv][ensgid] = self.typecast(
                                strip_value(arr[colidx]), self.field_type[colidx])

    def load_range_data(self):
        for line in file_util.gzopen(self.fname):
            line = line.decode('UTF-8')
            if line[0] != '#':
                arr = line.split('\t')
                arr[-1] = arr[-1].strip()
                chrom = arr[0].strip()
                spos = int(arr[1].strip())
                epos = int(arr[2].strip())
                d = {'spos': spos, 'epos': epos, 'v': []}
                for cidx in self.target_colidx:
                    if cidx == -999:
                        d['v'].append('')
                    else:
                        d['v'].append(strip_value(arr[cidx]))
                    try:
                        self.data[chrom].append(d)
                    except KeyError:
                        self.data[chrom] = [d]

    def get_range_data(self, geneinfo):
        rst = []
        geneinfo['epos'] = int(geneinfo['epos'])
        geneinfo['spos'] = int(geneinfo['spos'])
        try:
            for v1 in self.data[geneinfo['chrom']]:
                if geneinfo['epos'] >= v1['spos'] and geneinfo['spos'] <= v1['epos']:
                    rst.append(v1['v'])
                elif geneinfo['epos'] < v1['spos']:
                    break
        except KeyError:
            pass
        return rst

    def is_empty_value(self, v1):
        flag = False
        if v1 == "":
            flag = True
        if type(v1) == list:
            if len(v1) == 0:
                flag = True
            if len(v1) == 1 and v1[0] == "":
                flag = True
        return flag

    def get_data_dict_ensgid(self, ensgid, geneinfo, colname_type='name'):
        rstdata = []
        if self.is_range_data:
            datalist = self.get_range_data(geneinfo)
        else:
            try:
                datalist = self.data[ensgid]
            except KeyError:
                datalist = []

        if self.is_range_data and len(self.target_colnames) == 1:
            # for special case (CYTOBAND)
            d = []
            for datarow in datalist:
                d.append(datarow[0])
            if not self.is_empty_value(d):
                rstdata.append({self.target_colnames[0]: d})
        else:
            for datarow in datalist:
                d = {}
                for k in range(len(self.target_colnames)):
                    if colname_type == 'name2':
                        k1 = self.target_colnames2[k]
                    else:
                        k1 = self.target_colnames[k]
                    if not self.is_empty_value(datarow[k]) and datarow[k] is not None:
                        d[k1] = datarow[k]
                rstdata.append(d)

        return rstdata


class GeneSourceMerger():
    def __init__(self, source_readers):
        self.source_readers = source_readers
        self.ensgid_list = []
        self.colnames = []
        self.colnames2 = []
        self.set_colnames()
        self.merge_ensgid_list()

    def merge_ensgid_list(self):
        ensgid_map = {}
        for sid in self.source_readers.keys():
            s1 = self.source_readers[sid]
            # print(len(s1.ensgid_list), s1.ensgid_list[:3])
            for ensgid in s1.ensgid_list:
                if ensgid != '':
                    ensgid_map[ensgid] = 1
        self.ensgid_list = list(ensgid_map.keys())

    def set_colnames(self):
        self.colnames = {}
        self.colnames2 = {}
        for sid in self.source_readers.keys():
            s1 = self.source_readers[sid]
            self.colnames[sid] = s1.target_colnames
            self.colnames2[sid] = s1.target_colnames2

    def get_data_dict(self, ensgid):
        data = {}
        resv = self.get_reserved_data(ensgid)

        data = {}
        if 'spos' in resv.keys() and 'epos' in resv.keys():
            for r1 in RESERVED_COL:
                try:
                    v1 = resv[r1]
                except KeyError:
                    v1 = ''
                data[r1] = v1

            for sid in self.source_readers.keys():
                s1 = self.source_readers[sid]

                d1 = s1.get_data_dict_ensgid(ensgid, resv, 'name2')
                if s1.subembedded == "":
                    # if len(d1) > 2:
                    #     print(d1)
                    if len(d1) > 0:
                        data.update(d1[0])
                else:
                    if len(d1) > 0:
                        data.update({s1.subembedded: d1})
        return data

    def get_reserved_data(self, ensgid):
        resv = {}
        for sid in self.source_readers.keys():
            s1 = self.source_readers[sid]
            if s1.is_ensemblegene:
                for r1 in s1.reserved_data.keys():
                    try:
                        resv[r1] = s1.reserved_data[r1][ensgid]
                    except KeyError:
                        pass
        return resv

    def get_data_list(self, ensgid):
        data = []
        resv = self.get_reserved_data(ensgid)

        if 'spos' in resv.keys() and 'epos' in resv.keys():
            for r1 in RESERVED_COL:
                try:
                    r2 = resv[r1]
                except KeyError:
                    r2 = ''
                data.append(r2)

            for sid in self.source_readers.keys():
                s1 = self.source_readers[sid]
                if s1.is_range_data:
                    rdata = self.get_range_data(s1.data, resv)
                    if len(rdata) > 0:
                        data.extend(rdata[0])
                    else:
                        for cidx in s1.target_colidx:
                            data.append('')
                else:
                    try:
                        data.extend(s1.data[ensgid])
                    except KeyError:
                        for cidx in s1.target_colidx:
                            data.append('')
        print('get_data_list', data)
        return data


class GeneDataSourceFile():
    def __init__(self, opt):
        self.vartype = opt.vartype
        self.outtype = opt.outtype
        self.set_outfile_extension(opt.out, opt.outtype)
        self.datastruct = file_util.load_json(opt.ds)
        self.region = opt.region
        self.blocksize = opt.blocksize

    def set_outfile_extension(self, out, outtype):
        out2 = out.replace('.'+outtype+'.gz', '.' + outtype)
        if out2[-4:] != "." + outtype:
            out2 += "." + outtype
        self.out = out2

    def get_bed_header(self):
        header = []
        for r1 in RESERVED_COL:
            header.append(r1)
        info_arr = []
        for s1 in self.datastruct['gene_source']:
            if is_available(s1):
                subembedded = struct_util.get_dict_value(s1, "subembedded", "")

                info_header = []
                for f1 in s1['fields']:
                    if is_available(f1):
                        if 'name2' in f1.keys() and f1['name2'] != '':
                            fname = f1['name2']
                        else:
                            fname = f1['name']
                        if fname not in header:
                            info_header.append(s1['name'] + "_" + fname)
                info_arr.append('|'.join(info_header))
        header.append('\t'.join(info_arr))
        return '#' + '\t'.join(header)

    def save_json_out(self, data):
        # with open(self.out, 'w') as fp:
        fp = open(self.out, 'w')
        json.dump(data, fp)

    def save_file_header(self, fp):
        if self.outtype == "json":
            fp.write('[')

    def save_file_tail(self, fp):
        if self.outtype == "json":
            fp.write(']')

    def save_file_record(self, fp, data, prev_flag):

        if self.outtype == "json":
            flag = True
            if "CODING_GENE" in self.vartype:
                if data['gene_biotype'] not in ["protein_coding", "miRNA", "polymorphic_pseudogene"]:
                    flag = False
            if "MAIN_CHROM" in self.vartype:
                if data['chrom'] not in seq_util.MAIN_CHROM_LIST:
                    flag = False

            if flag:
                if prev_flag:
                    fp.write(',')
                # fp.write(json.dumps(data))
                fp.write(json.dumps(data, indent=2))
        else:
            fp.write(data)
        return flag

    # init function
    def make_single_source_file(self):
        file_util.check_dir(self.out)
        fp = open(self.out, 'w')
        # fp.write(self.get_bed_header() + '\n')

        sidlist = []
        source_readers = {}
        for s1 in self.datastruct['gene_source']:
            if is_available(s1):
                sid = s1['name']
                sidlist.append(sid)

                datafile_path = ''
                if 'datafile_path' in self.datastruct.keys():
                    datafile_path = self.datastruct['datafile_path']
                source_readers[sid] = GeneSourceReader(s1, datafile_path)

        gsm = GeneSourceMerger(source_readers)

        i = 0
        self.save_file_header(fp)
        flag_write = False
        for ensgid in gsm.ensgid_list:
            # data = gsm.get_data_list(ensgid)
            # line = '\t'.join(data)
            data = gsm.get_data_dict(ensgid)
            if len(data.keys()) > 0:

                flag2 = self.save_file_record(fp, data, flag_write)
                if flag2:
                    flag_write = True
                i += 1

            if i > 1000:
                pass
                # break
        self.save_file_tail(fp)
        print('Saved', self.out)
        fp.close()

        # file check
        print('Json checking...', self.out)
        fp = open(self.out)
        data = json.load(fp)
