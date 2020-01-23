import file_util
import vcf_util
import tabix
import time


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


class TSVBlockReader():

    def __init__(self, fname):
        self.fname = fname
        self.target_colidx = []
        self.header = []
        self.filter_start_with = {}
        self.filter_skip_equal_str = {}
        self.delimiter = {}
        self.tp = None
        self.chrom = ''
        self.read_header()

    def read_header(self):
        for line in file_util.gzopen(self.fname.replace('#CHROM#', "1")):
            line = line.decode('UTF-8')
            if line[0] == '#' and line[1] != '#':
                self.header = line[1:].strip().split('\t')
                break

    def set_target_column(self, fields_structure):
        for f1 in fields_structure:
            if f1['name'] in self.header:
                self.target_colidx.append(self.header.index(f1['name']))
                if 'start_with' in f1.keys():
                    self.filter_start_with[self.header.index(f1['name'])] = f1['start_with']
                if 'skip_equal_str' in f1.keys():
                    self.filter_skip_equal_str[self.header.index(f1['name'])] = f1['skip_equal_str']
                if 'delimiter' in f1.keys():
                    self.delimiter[self.header.index(f1['name'])] = f1['delimiter']

    def add_block(self, block, regionstr, sid):
        region = pars_region_str(regionstr)
        if self.tp is None:
            # FIXME: temporary code
            # if 'VEP' in self.fname:
            #     fname = "/n/data1/hms/dbmi/park/daniel/BiO/Research/mutanno/PRECALVEP/chr" + region['chrom']
            #     fname += "/" + str(int(region['spos'] / 1000000))
            #     fname += "/chr" + self.ori_region.replace(':', '_').replace('-', '_') + ".tsv.gz"
            #     print(fname)
            #     self.fname = fname
            sourcefile = self.fname.replace('#CHROM#', region['chrom'])
            if file_util.is_exist(sourcefile) and file_util.is_exist(sourcefile + ".tbi"):
                self.tp = tabix.open(self.fname.replace('#CHROM#', region['chrom']))
        if self.tp is not None:
            for arr in self.tp.querys(regionstr):
                pos = int(arr[1])
                if pos >= region['spos'] and pos <= region['epos']:
                    refalt = arr[self.header.index('REF')] + '_' + arr[self.header.index('ALT')]
                    filter = True
                    cont = ''
                    for cidx in self.target_colidx:
                        if cont != '':
                            cont += '|'
                        delimiter = ''
                        if cidx in self.delimiter.keys():
                            delimiter = self.delimiter[cidx]
                        if cidx >= 0:
                            cont += vcf_util.encode_infovalue(arr[cidx], delimiter)

                        if cidx in self.filter_start_with.keys():
                            filter = arr[cidx].startswith(self.filter_start_with[cidx])
                        if cidx in self.filter_skip_equal_str.keys():
                            filter = not (arr[cidx] == self.filter_skip_equal_str[cidx])
                    if filter:
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
            # FIXME: temporary code
            # self.block_readers[sid].ori_region = self.ori_region
            self.block = self.block_readers[sid].add_block(self.block, range, sid)

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
            # for chrom in seq_util.MAIN_CHROM_LIST:
            sid = s1['name']
            sidlist.append(sid)
            blockreaders[sid] = TSVBlockReader(s1['datafile'])
            blockreaders[sid].set_target_column(s1['fields'])

        tbm = TSVBlockMerger(blockreaders, self.region, self.blocksize)

        while not tbm.eof:
            start = time.time()
            block = tbm.get_block_tsi()
            fp.write(block)
            lastpos = str(tbm.block_lpos)
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
