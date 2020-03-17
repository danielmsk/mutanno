import os
import tabix
import time
import file_util
import vcf_util
import struct_util
import external_functions

REFIDX = 3
ALTIDX = 4
entrezmap = {}
refseqmap = {}


def pars_region_str(region):
    r = None
    if region != "":
        arr = region.split(':')
        arr2 = arr[1].split('-')
        r = {}
        r['chrom'] = arr[0].replace('chr', '')
        r['spos'] = int(arr2[0])
        r['epos'] = int(arr2[1])
    return r


def is_available(field):
    return struct_util.is_available(field)


class TSVBlockReader():

    def __init__(self, datastruct, datafile_path):
        self.datastruct = datastruct
        self.fname = os.path.join(datafile_path, self.datastruct['datafile'])
        self.fileformat = self.datastruct['format']
        self.target_colidx = []
        self.header = []
        self.filter_start_with = {}
        self.filter_skip_equal_str = {}
        self.delimiter = {}
        self.defaultvalue = {}
        self.field_function = {}
        self.field_function_param = {}
        self.field_names = []
        self.tp = None
        self.chrom = ''
        self.chrompre = ''
        self.altidx = ALTIDX
        self.refidx = REFIDX

        if 'ref_column_index' in self.datastruct.keys():
            self.refidx = self.datastruct['ref_column_index']
        if 'ref_column_index' in self.datastruct.keys():
            self.altidx = self.datastruct['alt_column_index']
        if 'chrompre' in self.datastruct.keys():
            self.chrompre = self.datastruct['chrompre']

        self.read_header()
        self.set_target_column(self.datastruct['fields'])

    def read_header(self):
        for line in file_util.gzopen(self.fname.replace('#CHROM#', "1")):
            line = line.decode('UTF-8')
            if line[0] == '#' and line[1] != '#':
                if self.fileformat == "tsv" or self.fileformat == "bed" or self.fileformat == "vcf":
                    self.header = line[1:].strip().split('\t')
                if self.fileformat == "tsi":
                    self.header = line[1:].strip().split('\t')[-1].split('|')
                break

    def set_target_column(self, fields_structure):
        cidx_no = 0
        for f1 in fields_structure:
            if is_available(f1):
                self.field_names.append(f1['name'])
                if f1['name'] in self.header:
                    self.target_colidx.append(self.header.index(f1['name']))
                    if 'start_with' in f1.keys():
                        self.filter_start_with[self.header.index(f1['name'])] = f1['start_with']
                    if 'skip_equal_str' in f1.keys():
                        self.filter_skip_equal_str[self.header.index(f1['name'])] = f1['skip_equal_str']
                    if 'delimiter' in f1.keys():
                        self.delimiter[self.header.index(f1['name'])] = f1['delimiter']
                    if 'default' in f1.keys():
                        self.defaultvalue[self.header.index(f1['name'])] = f1['default']
                else:
                    self.target_colidx.append(-999)

                cidx_no += 1
                if 'function' in f1.keys() and f1['function'] != '':
                    self.field_function[cidx_no] = f1['function']
                    self.field_function_param[cidx_no] = ''
                    if 'param' in f1.keys() and f1['param'] != '':
                        self.field_function_param[cidx_no] = f1['param']

    def add_block_bed(self, block, regionstr, sid):
        region = pars_region_str(regionstr)
        # print(sid, self.fname, region)
        if self.tp is None:
            sourcefile = self.fname.replace('#CHROM#', region['chrom'])
            if file_util.is_exist(sourcefile) and file_util.is_exist(sourcefile + ".tbi"):
                self.tp = tabix.open(self.fname.replace('#CHROM#', region['chrom']))
        if self.tp is not None:
            bedblock = []
            for arr in self.tp.querys(regionstr):
                d = {}
                d['spos'] = int(arr[1])
                d['epos'] = int(arr[2])
                for cidx in self.target_colidx:
                    if cidx == -999:
                        d[cidx] = ''
                    else:
                        d[cidx] = arr[cidx]
                bedblock.append(d)

            if len(bedblock) > 0:
                for pos in block.keys():
                    cont = ''
                    for d in bedblock:
                        if pos > d['spos'] and pos <= d['epos']:
                            if cont != '':
                                cont += ','
                            variantkey = region['chrom'] + '_' + str(pos)
                            flag_filter, cont2 = self.get_selected_fields_in_block(d, variantkey)
                            cont += cont2
                    if cont != '':
                        for refalt in block[pos].keys():
                            block[pos][refalt][sid] = cont
        return block

    def get_selected_fields_in_block(self, arr, variantkey):
        cidx_no = 0
        flag_filter = True
        arr_selected_fields = []
        for cidx in self.target_colidx:
            cidx_no += 1
            if cidx == -999:
                cidxvalue = ''
            else:
                cidxvalue = arr[cidx]

            if cidx_no in self.field_function.keys():
                # print(arr_selected_fields)
                exec_str = "external_functions." + self.field_function[cidx_no]
                exec_str += '('
                paramstr_list = []
                for param in self.field_function_param[cidx_no].split(','):
                    if param == "mutanno_value_variantkey":
                        pstr = 'variantkey'
                    elif cidx_no-1 == self.field_names.index(param):
                        pstr = 'cidxvalue'
                    else:
                        pstr = 'arr_selected_fields['
                        pstr += str(self.field_names.index(param))
                        pstr += ']'
                    paramstr_list.append(pstr)
                exec_str += ','.join(paramstr_list)
                exec_str += ')'
                # print(exec_str)
                cidxvalue = eval(exec_str)

            delimiter = ''
            if cidx in self.delimiter.keys():
                delimiter = self.delimiter[cidx]

            arr_selected_fields.append(vcf_util.encode_infovalue(cidxvalue, delimiter))

            if cidx in self.filter_start_with.keys():
                flag_filter = cidxvalue.startswith(self.filter_start_with[cidx])
            if cidx in self.filter_skip_equal_str.keys():
                flag_filter = not (cidxvalue == self.filter_skip_equal_str[cidx])

        cont = '|'.join(arr_selected_fields)
        return flag_filter, cont

    def add_block(self, block, regionstr, sid):
        region = pars_region_str(regionstr)
        if self.tp is None:
            sourcefile = self.fname.replace('#CHROM#', region['chrom'])
            if file_util.is_exist(sourcefile) and file_util.is_exist(sourcefile + ".tbi"):
                self.tp = tabix.open(self.fname.replace('#CHROM#', region['chrom']))
        if self.tp is not None:
            # print(self.chrompre + regionstr)
            try:
                # print(attributes(self.tp))
                # if True:
                for arr in self.tp.querys(self.chrompre + regionstr):
                    pos = int(arr[1])
                    if pos >= region['spos'] and pos <= region['epos']:
                        refalt = arr[self.refidx] + '_' + arr[self.altidx]

                        if self.fileformat == 'tsi':
                            if "," in arr[-1]:
                                arrsection = []
                                for f1 in arr[-1].split(','):
                                    arrsection.append(f1.split('|'))
                            else:
                                arrsection = [arr[-1].split('|')]
                        else:
                            arrsection = [arr]

                        for sec in arrsection:
                            variantkey = region['chrom'] + '_' + \
                                str(pos) + '_' + arr[self.refidx] + '_' + arr[self.altidx]
                            flag_filter, cont = self.get_selected_fields_in_block(sec, variantkey)
                            if flag_filter:
                                try:
                                    block[pos]
                                except KeyError:
                                    block[pos] = {}
                                try:
                                    block[pos][refalt]
                                except KeyError:
                                    block[pos][refalt] = {}
                                try:
                                    block[pos][refalt][sid] += ',' + cont
                                except KeyError:
                                    block[pos][refalt][sid] = cont
            except tabix.TabixError:
                # there is no annotation in the chrom. (ex. Y, M)
                pass
        return block


class TSVBlockMerger():
    def __init__(self, block_readers, region, block_size):
        self.block_readers = block_readers
        self.block_size = block_size
        self.region = pars_region_str(region)
        self.ori_region = region
        self.block_spos = 0
        self.block_epos = 0
        self.block_lpos = 0
        self.block = {}
        self.eof = False

    def read_block(self):
        if self.block_spos == 0:
            self.block_spos = self.region['spos']
        else:
            self.block_spos = self.block_epos + 1

        if (self.block_spos + self.block_size) <= self.region['epos']:
            self.block_epos = self.block_spos + self.block_size - 1
        else:
            self.block_epos = self.region['epos']
        range = self.region['chrom'] + ':' + str(self.block_spos) + '-' + str(self.block_epos)

        self.block = {}
        for sid in self.block_readers.keys():
            if self.block_readers[sid].fileformat != "bed":
                self.block = self.block_readers[sid].add_block(self.block, range, sid)
        for sid in self.block_readers.keys():
            if self.block_readers[sid].fileformat == "bed":
                self.block = self.block_readers[sid].add_block_bed(self.block, range, sid)

    def get_block_tsi(self):
        self.read_block()
        if self.block_epos >= self.region['epos']:
            self.eof = True

        cont = ''
        poslist = list(self.block.keys())
        for pos in sorted(poslist):
            for refalt in self.block[pos].keys():
                info = ''
                for sid in self.block[pos][refalt].keys():
                    if info != '':
                        info += ';'
                    info += sid + '=' + self.block[pos][refalt][sid]
                cont += self.region['chrom'] + '\t' + str(pos) + '\t\t' + refalt.replace('_', '\t') + '\t' + info + '\n'
            self.block_lpos = pos
        return cont


class DataSourceFile():
    def __init__(self, opt):
        self.set_outfile_extension(opt['out'])
        self.datastruct = file_util.load_json(opt['ds'])
        self.region = opt['region']
        self.blocksize = opt['blocksize']

    def set_outfile_extension(self, out):
        out2 = out.replace('.tsi.gz', '.tsi')
        if out2[-4:] != ".tsi":
            out2 += ".tsi"
        self.out = out2

    def get_tsv_header(self):
        header = ["#CHROM", "POS", "ID", "REF", "ALT"]
        for s1 in self.datastruct['source']:
            for f1 in s1['fields']:
                if is_available(f1):
                    if 'name2' in f1.keys() and f1['name2'] != '':
                        header.append(f1['name2'])
                    else:
                        header.append(f1['name'])
        return '\t'.join(header)

    def get_tsi_header(self):
        header = ["#CHROM", "POS", "ID", "REF", "ALT"]
        info_arr = []
        for s1 in self.datastruct['source']:
            info_header = []
            for f1 in s1['fields']:
                if is_available(f1):
                    if 'name2' in f1.keys() and f1['name2'] != '':
                        info_header.append(f1['name2'])
                    else:
                        info_header.append(f1['name'])
            info_arr.append(s1['name'] + '=' + '|'.join(info_header))
        header.append(';'.join(info_arr))
        return '\t'.join(header)

    def make_single_source_file(self):
        file_util.check_dir(self.out)
        fp = open(self.out, 'w')
        fp.write(self.get_tsi_header() + '\n')

        sidlist = []
        blockreaders = {}
        for s1 in self.datastruct['source']:
            sid = s1['name']
            sidlist.append(sid)
            datafile_path = ''
            if 'datafile_path' in self.datastruct.keys():
                datafile_path = self.datastruct['datafile_path']
            blockreaders[sid] = TSVBlockReader(s1, datafile_path)

        tbm = TSVBlockMerger(blockreaders, self.region, self.blocksize)

        while not tbm.eof:
            start = time.time()
            block = tbm.get_block_tsi()
            fp.write(block)
            # lastpos = str(tbm.block_lpos)
            end = time.time()
            elapsed = end - start
            log = 'processing.. ' + str(round(elapsed, 3)) + ' sec elapsed. '
            log += str(len(tbm.block.keys())) + ' variants added. '
            log += str(tbm.block_spos) + '~' + str(tbm.block_epos)
            print(log)
        print('Saved', self.out)
        fp.close()

        # FIXME: convert using bgzip lib
        # time.sleep(3)
        # proc_util.run_cmd('tabixgz ' + self.out)
        # time.sleep(3)
        # for line in file_util.gzopen(self.out + '.gz'):
        #     line = line.decode('UTF-8')
        #     lastposgz = line.split('\t')[1]
        #
        # if (lastpos == lastposgz and file_util.is_exist(self.out + '.gz.tbi')) or lastposgz == "POS":
        #     cmd = "rm -rf " + self.out + ";"
        #     cmd += "rm -rf " + self.out + ".gz.tbi;"
        #     proc_util.run_cmd(cmd)
