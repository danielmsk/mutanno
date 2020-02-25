import os
import time
import file_util
import struct_util


RESERVED_COL = ["chrom", "spos", "epos", "ensgid"]

def is_available(field):
    return struct_util.is_available(field)

def is_reserved_column(f1):
    flag = False
    if f1['name'] in RESERVED_COL or ('name2' in f1.keys() and f1['name2'] in RESERVED_COL):
        flag = True
    return flag

def strip_value(v1):
    v1 = v1.replace('"','')
    return v1

class GeneSourceReader():
    def __init__(self, datastruct, datafile_path):
        self.datastruct = datastruct
        self.source_name = self.datastruct['name']
        self.fname = os.path.join(datafile_path, self.datastruct['datafile'])
        self.fileformat = self.datastruct['format']
        self.target_colidx = []
        self.header = []
        self.field_function = {}
        self.field_function_param = {}
        self.field_names = []
        self.key_colidx = -999
        self.reserved_colidx = {}
        self.reserved_data = {}
        self.data = {}
        self.ensgid_list = []
        self.is_range_data = False
        self.is_ensemblegene = struct_util.get_dict_value(datastruct,'is_ensemblegene',False)
        self.read_header()
        print("loading..", self.source_name)
        # print(self.header)
        self.set_target_column(self.datastruct['fields'])

        if self.key_colidx == -999:
            self.load_range_data()
            self.is_range_data = True
        else:
            self.load_data()

    def read_header(self):
        for line in file_util.gzopen(self.fname):
            line = line.decode('UTF-8')
            if line[0] == '#' and line[1] != '#':
                self.header = line[1:].strip().split('\t')
                break

    def set_target_column(self, fields_structure):
        cidx_no = 0
        for f1 in fields_structure:

            if is_available(f1):
                self.field_names.append(f1['name'])

                # print("self.header",self.header, f1['name'], self.key_colidx)
                if f1['name'] in self.header:
                    if is_reserved_column(f1):
                        for resv in RESERVED_COL:
                            if f1['name'] == resv or ('name2' in f1.keys() and f1['name2'] == resv):
                                self.reserved_colidx[resv] = self.header.index(f1['name'])
                    else:
                        self.target_colidx.append(self.header.index(f1['name']))
                else:
                    if not is_reserved_column(f1):
                        self.target_colidx.append(-999)

                cidx_no += 1
                if 'function' in f1.keys() and f1['function'] != '':
                    self.field_function[cidx_no] = f1['function']
                    self.field_function_param[cidx_no] = ''
                    if 'param' in f1.keys() and f1['param'] != '':
                        self.field_function_param[cidx_no] = f1['param']

        if 'ensgid' in self.reserved_colidx.keys():
            self.key_colidx = self.reserved_colidx['ensgid']


    def rm_version_in_ensgid(self, ensg):
        return ensg.strip().split('.')[0]

    def load_data(self):
        for line in file_util.gzopen(self.fname):
            line = line.decode('UTF-8')
            if line[0] != '#':
                arr = line.split('\t')
                arr[-1] = arr[-1].strip()
                ensgid = self.rm_version_in_ensgid(arr[self.key_colidx])
                self.ensgid_list.append(ensgid)
                self.data[ensgid] = []
                for cidx in self.target_colidx:
                    if cidx == -999:
                        v1 = ""
                    else:
                        v1 = strip_value(arr[cidx])
                    self.data[ensgid].append(v1)

                for resv in self.reserved_colidx.keys():
                    try:
                        tmp = self.reserved_data[resv]
                    except KeyError:
                        self.reserved_data[resv] = {}

                    try:
                        tmp = self.reserved_data[resv][ensgid]
                    except KeyError:
                        if resv == 'ensgid':
                            self.reserved_data[resv][ensgid] = self.rm_version_in_ensgid(arr[self.reserved_colidx[resv]])
                        else:
                            self.reserved_data[resv][ensgid] = strip_value(arr[self.reserved_colidx[resv]])

    def load_range_data(self):
        for line in file_util.gzopen(self.fname):
            line = line.decode('UTF-8')
            if line[0] != '#':
                arr = line.split('\t')
                arr[-1] = arr[-1].strip()
                chrom = arr[0].strip()
                spos = int(arr[1].strip())
                epos = int(arr[2].strip())
                d = {'spos':spos, 'epos':epos, 'v':[]}
                for cidx in self.target_colidx:
                    if cidx == -999:
                        d['v'].append('')
                    else:
                        d['v'].append(strip_value(arr[cidx]))
                    try:
                        self.data[chrom].append(d)
                    except KeyError:
                        self.data[chrom] = [d]

class GeneSourceMerger():
    def __init__(self, source_readers):
        self.source_readers = source_readers
        self.ensgid_list = []
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

    def get_range_data(self, datalist, geneinfo):
        rst = []
        geneinfo['epos'] = int(geneinfo['epos'])
        geneinfo['spos'] = int(geneinfo['spos'])
        try:
            for v1 in datalist[geneinfo['chrom']]:
                if geneinfo['epos'] >= v1['spos'] and geneinfo['spos'] <= v1['epos']:
                    rst.append(v1['v'])
                elif geneinfo['epos'] < v1['spos']:
                    break
        except KeyError:
            pass
        return rst

    def get_line(self, ensgid):
        cont = []
        rst = ""
        resv = {}
        for sid in self.source_readers.keys():
            s1 = self.source_readers[sid]
            if s1.is_ensemblegene:
                for r1 in s1.reserved_data.keys():
                    try:
                        resv[r1] = s1.reserved_data[r1][ensgid]
                    except KeyError:
                        pass

        if 'spos' in resv.keys() and 'epos' in resv.keys():
            for r1 in RESERVED_COL:
                try:
                    r2 = resv[r1]
                except KeyError:
                    r2 = ''
                cont.append(r2)

            for sid in self.source_readers.keys():
                s1 = self.source_readers[sid]
                if s1.is_range_data:
                    rdata = self.get_range_data(s1.data, resv)
                    if len(rdata) > 0:
                        cont.extend(rdata[0])
                    else:
                        for cidx in s1.target_colidx:
                            cont.append('')
                else:
                    try:
                        cont.extend(s1.data[ensgid])
                    except KeyError:
                        for cidx in s1.target_colidx:
                            cont.append('')
            rst = '\t'.join(cont) + '\n'
        return rst


class GeneDataSourceFile():
    def __init__(self, opt):
        self.set_outfile_extension(opt['out'])
        self.datastruct = file_util.load_json(opt['ds'])
        self.region = opt['region']
        self.blocksize = opt['blocksize']

    def set_outfile_extension(self, out):
        out2 = out.replace('.bed.gz', '.bed')
        if out2[-4:] != ".bed":
            out2 += ".bed"
        self.out = out2

    def get_bed_header(self):
        header = [] 
        for r1 in RESERVED_COL:
            header.append(r1)
        info_arr = []
        for s1 in self.datastruct['gene_source']:
            if is_available(s1):
                info_header = []
                for f1 in s1['fields']:
                    if is_available(f1):
                        if 'name2' in f1.keys() and f1['name2'] != '':
                            fname = f1['name2']
                        else:
                            fname = f1['name']
                        if fname not in header:
                            info_header.append(s1['name'] + "_" + fname)
                info_arr.append('\t'.join(info_header))
        header.append('\t'.join(info_arr))
        return '#' + '\t'.join(header)

    # init function
    def make_single_source_file(self):
        file_util.check_dir(self.out)
        fp = open(self.out, 'w')
        fp.write(self.get_bed_header() + '\n')

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
        for ensgid in gsm.ensgid_list:
            line = gsm.get_line(ensgid)
            if line != "":
                fp.write(line)
            i += 1
            if i > 1000:
                pass
                # break
            
        print('Saved', self.out)
        fp.close()




