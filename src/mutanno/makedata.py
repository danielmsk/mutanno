import os
import tabix
import time
from .util import file_util
from .util import vcf_util
from .util import proc_util
from .util import struct_util
from .util import region_util
from .model import datasource
from .annot2.annotvcf import VCFAnnotator, VCFVariant
from .model.datastructure import DataSourceListStructure
from .model.datasource import DataSourceList
from . import external_functions
from .renderer import TSIRenderer

REFIDX = 3
ALTIDX = 4
entrezmap = {}
refseqmap = {}



class BlankVCFReader():
    def __init__(self, opt, dslist):
        self.dslist = dslist
        self.opt = opt
        self.region = region_util.Region(opt.region)
        self.tpos = self.region.spos
        self.variantkeylist = []
        self.eof = False
        self.colnames = vcf_util.VCF_COL

    def get_variant(self):
        if len(self.variantkeylist) == 0 and self.tpos < self.region.epos:
            self.load_variantkeys()

        if len(self.variantkeylist) > 0:
            vkey = self.variantkeylist.pop(0)
            vcfline = vkey.replace('_', '\t') + '\t\t\t\t'
            variant = VCFVariant(vcfline, self.colnames, self.opt)
            self.tpos = int(vkey.split('_')[1])
        else:
            variant = None
            self.eof = True
        return variant

    def load_variantkeys(self):
        for r1 in self.region.get_split_region(1000):
            variantkeys = self.dslist.get_variantkeys(r1.chrom, r1.spos, r1.epos)
            poslist = list(variantkeys.keys())
            for pos in sorted(poslist):
                for vkey in variantkeys[pos]:
                    self.variantkeylist.append(vkey)


class DataSourceGenerator(VCFAnnotator):
    def __init__(self, opt):
        self.opt = opt
        self.dslist = DataSourceList(self.opt.ds)
        self.region = region_util.Region(opt.region)
        self.out = None
        self.set_outfile_extension()
        self.vcfreader = BlankVCFReader(self.opt, self.dslist)
        self.tsirenderer = TSIRenderer()

    def set_outfile_extension(self):
        out2 = self.opt.out.replace('.tsi.gz', '.tsi')
        if out2[-4:] != ".tsi":
            out2 += ".tsi"
        self.out = out2

    def get_tsi_header(self):
        header = ["#CHROM", "POS", "ID", "REF", "ALT"]
        info_arr = []
        for s1 in self.dslist.available_source_list:
            info_header = []
            for f1 in s1.available_field_list:
                info_header.append(f1.name2)
            info_arr.append(s1.name + '=' + '|'.join(info_header))
        header.append(';'.join(info_arr))
        return '\t'.join(header)

    def print_i(self, ivar, variant, stime1, stime0):
        if ivar % 500 == 0:
            etime1 = time.time()
            elapsed1 = etime1 - stime1
            elapsed0 = etime1 - stime0
            log = 'processed... ' + str(ivar) 
            log += " elapsed:" + str(round(elapsed1, 2)) + "s " + str(round(elapsed0, 2)) + "s , time for 9G variants:" + str(round(elapsed0 / ivar * 9000000000 / 3600, 3)) + "hr"
            # print(log, ivar, variant)

    def make_single_source_file(self):
        stime0 = time.time()
        self.open_outpointer()
        self.fp.write(self.get_tsi_header() + '\n')

        ivar = 0
        stime1 = time.time()
        while(not self.vcfreader.eof):
            ivar += 1
            variant = self.vcfreader.get_variant()
            if self.vcfreader.eof or variant == None:
                break
            variant.set_annot(self.dslist)
            self.fp.write(self.tsirenderer.render_vcfvariant(variant, False, True))
            self.print_i(ivar, variant, stime1, stime0)
            stime1 = time.time()
        
        self.close_outpointer()
        etime0 = time.time()
        elapsed0 = etime0 - stime0
        print("total:", elapsed0, ", time per variant:", elapsed0 / ivar, ", time for 9G variants:", str(round(elapsed0 / ivar * 9000000000 / 3600, 3)) + "hr")
        # time.sleep(10)
        # proc_util.run_cmd('gz ' + self.out)
        

    def make_single_source_file_with_block(self):
        fp = open(self.opt.out, 'w')
        blockreaders = {}
        for s1 in self.dslist.available_source_list:
            blockreaders[s1.name] = TSVBlockReader(s1)

        tbm = TSVBlockMerger(blockreaders, self.region, int(self.opt.blocksize))
        print_header = True
        while not tbm.eof:
            start = time.time()
            block = tbm.get_block_tsi()
            if print_header:
                fp.write(tbm.get_header())
                print_header = False
            fp.write(block)
            
            end = time.time()
            elapsed = end - start
            log = 'processing.. ' + str(round(elapsed, 3)) + ' sec elapsed. '
            log += str(len(tbm.block.keys())) + ' variants added. '
            log += str(tbm.block_spos) + '~' + str(tbm.block_epos)
            self.log(log)
        fp.close()
        self.log('Saved.. ' + self.out)

    def log(self, msg):
        file_util.fileSave(self.out + '.log', msg + '\n', 'a')
        print(msg)

    def check_datasourcefile(self):
        mt_seq_reader = datasource.MutannoDataSequentialReader(self.out)

        dslist = datasource.DataSourceList()
        dslist.set_datastructure(self.datastruct)

        for mtdata in mt_seq_reader:
            annot = dslist.get_sourcedata_variant(mtdata['CHROM'], mtdata['POS'], mtdata['REF'], mtdata['ALT'])
            
            break
        pass



class TSVBlockReader():

    def __init__(self, datasource):
        self.source = datasource
        self.fileformat = self.source.format
        self.target_colidx = []
        self.header = []
        self.filter_start_with = {}
        self.filter_skip_equal_str = {}
        self.filter_keep_equal_str = {}
        self.delimiter = {}
        self.defaultvalue = {}
        self.field_function = {}
        self.field_function_param = {}
        self.field_names = []
        self.tp = None
        self.chrom = ''
        self.refidx = self.source.ref_column_index
        self.altidx = self.source.alt_column_index
 
    def add_block_bed(self, block, region, sid):
        if self.tp is None:
            sourcefile = self.source.sourcefile2.replace('#CHROM#', region.chrom)
            if file_util.is_exist(sourcefile) and file_util.is_exist(sourcefile + ".tbi"):
                self.tp = tabix.open(sourcefile)
                self.source.set_header(sourcefile)

        if self.tp is not None:
            self.set_target_column()
            for arr in self.tp.querys(region.region_str):
                d = {}
                spos = int(arr[1])
                epos = int(arr[2])
                ref = arr[self.refidx]
                if (spos+1) >= region.spos and (spos+1) <= region.epos:
                    cont = (spos,epos,ref, self.source.name + "=" + "|".join(self.select_field(arr, self.source.header)))
                    try:
                        block[spos]
                    except KeyError:
                        block[spos] = {}
                    try:
                        block[spos]['BED'].append(cont)
                    except KeyError:
                        block[spos]['BED'] = [cont]
        return block

    def set_target_column(self):
        for idx, h1 in enumerate(self.source.header):
            if self.source.format == "bed":
                self.target_colidx = list(range(3, len(self.source.header)))
            elif self.fileformat == "tsi":
                self.target_colidx = [len(self.source.header)-1]
            else:
                for idx, h1 in enumerate(self.source.header):
                    if h1.upper() not in ['CHROM','POS','REF','ALT','VEP']:
                        if idx not in self.target_colidx:
                            self.target_colidx.append(idx)

    def select_field(self, arr, header):
        rst = []
        for idx in self.target_colidx:
            rst.append(vcf_util.encode_value(arr[idx]))
        return rst


    def add_block(self, block, region, sid):
        if self.tp is None:
            sourcefile = self.source.sourcefile2.replace('#CHROM#', region.chrom)
            if file_util.is_exist(sourcefile) and file_util.is_exist(sourcefile + ".tbi"):
                self.tp = tabix.open(sourcefile)
                self.source.set_header(sourcefile)
        if self.tp is not None:
            self.set_target_column()
            try:
                for arr in self.tp.querys(self.source.chrompre + region.region_str):
                    pos = int(arr[1])
                    if pos >= region.spos and pos <= region.epos:
                        refalt = arr[self.refidx] + '_' + arr[self.altidx]
                        try:
                            block[pos]
                        except KeyError:
                            block[pos] = {}
                        try:
                            block[pos][refalt]
                        except KeyError:
                            block[pos][refalt] = {}

                        if self.fileformat == "tsi":
                            cont = arr[-1].replace('&', '~').replace(' ', '%20')
                        else:
                            cont = '|'.join(self.select_field(arr, self.source.header))

                        try:
                            block[pos][refalt][sid].append(cont) 
                        except KeyError:
                            block[pos][refalt][sid] = [cont]
            except tabix.TabixError:
                pass
        
        return block

    def get_header(self):
        rst = []
        for idx in self.target_colidx:
            rst.append(self.source.header[idx])
        return rst



class TSVBlockMerger():
    def __init__(self, block_readers, region, block_size):
        self.block_readers = block_readers
        self.block_size = block_size
        self.region = region  # region object
        self.block_spos = 0
        self.block_epos = 0
        self.block_lpos = 0
        self.block = {}
        self.eof = False

    def read_block(self):
        if self.block_spos == 0:
            self.block_spos = self.region.spos
        else:
            self.block_spos = self.block_epos + 1

        if (self.block_spos + self.block_size) <= self.region.epos:
            self.block_epos = self.block_spos + self.block_size - 1
        else:
            self.block_epos = self.region.epos
        this_range = region_util.Region(self.region.chrom + ':' + str(self.block_spos) + '-' + str(self.block_epos))


        self.block = {}
        for sid in self.block_readers.keys():
            if self.block_readers[sid].fileformat != "bed":
                self.block = self.block_readers[sid].add_block(self.block, this_range, sid)

        for sid in self.block_readers.keys():
            if self.block_readers[sid].fileformat == "bed":
                self.block = self.block_readers[sid].add_block_bed(self.block, this_range, sid)

    def get_block_tsi(self):
        self.read_block()
        if self.block_epos >= self.region.epos:
            self.eof = True

        cont = ''
        poslist = list(self.block.keys())
        for pos in sorted(poslist):
            for refalt in self.block[pos].keys():
                if refalt == "BED":
                    for bedline in self.block[pos]['BED']:
                        spos = bedline[0]
                        epos = bedline[1]
                        ref = bedline[2]
                        infoarr = bedline[3]
                        cont += self.region.chrom + '\t' + str(spos+1) + '\t\t' + ref + '\t\t\t\t' + 'END=' + str(epos) +  '\t' + infoarr + '\n'
                else:
                    infoarr = []
                    for sid in self.block_readers.keys():
                        if sid in self.block[pos][refalt].keys():
                            if self.block[pos][refalt][sid] != "":
                                infoarr.append(sid + '=' + ",".join(self.block[pos][refalt][sid]))
                                
                        cont += self.region.chrom + '\t' + str(pos) + '\t\t' + refalt.replace('_', '\t') + '\t\t\t\t' + ';'.join(infoarr) + '\n'
            

                
            self.block_lpos = pos
        return cont

    def get_header(self):
        cont = ""
        for sid in self.block_readers.keys():
            header = self.block_readers[sid].get_header()
            if cont != "":
                cont += ";"
            cont += sid + '=' + '|'.join(header)
        cont = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t" + cont+ "\n"
        return cont



