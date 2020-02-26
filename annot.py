#!/usr/bin/env python
# annot.py
# made by Daniel Minseok Kwon
#########################
import tabix
import file_util
import vcf_util
import time

VCFCOLIDX = {'CHROM': 0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4, 'QUAL': 5, 'FILTER': 6, 'INFO': 7, 'FORMAT': 8}
TSVCOLIDX = {'CHROM': 0, 'POS': 1, 'REF': 2, 'ALT': 3}
DATAHEADER = {}

MAXBUFF = 1000000


def load_entrez_refseq():
    path = "/home/mk446/bio/mutanno/DATASOURCE/ENSEMBL/hg38/"
    entrezmap = {}
    refseqmap = {}
    try:
        for line in file_util.gzopen(path + 'Homo_sapiens.GRCh38.98.entrez.tsv.gz'):
            line = line.decode('UTF-8')
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            entrezmap[arr[0].strip()] = arr[3].strip()
        
        for line in file_util.gzopen(path + 'Homo_sapiens.GRCh38.98.refseq.sorted.tsv.gz'):
            line = line.decode('UTF-8')
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            enst_id = arr[1].strip()
            refseqmap[enst_id] = arr[3].strip()
    except FileNotFoundError:
        pass
    return entrezmap, refseqmap


entrezmap, refseqmap = load_entrez_refseq()


def conv_infovalue(v1):
    if v1 == '.' or v1 == '-':
        v1 = ''
    v1 = v1.replace(";", "%3B").replace("=", "%3D").replace("|", "%7C").replace(",", "%2C")
    v1 = v1.replace('"', "%22").replace('~', "%7E").replace(' ', '%20')
    # v1 = urllib.parse.quote(v1)
    return v1


def get_dict_value(dict1, key1, default):
    rst = default
    if key1 in dict1.keys():
        rst = dict1[key1]
    return rst


def get_sourcename(s1):
    s1name2 = get_dict_value(s1, 'name2', '')
    s1name = get_dict_value(s1, 'name', '')
    if s1name2 != '':
        s1name = s1name2
    return s1name


class AnnotBlock():
    def __init__(self, source, tabixpointer, header, buff=1000):
        self.tabixpointer = tabixpointer
        self.source = source
        self.sourcename = source['name']
        self.header = header
        self.fileformat = source['format']
        self.chrompre = get_dict_value(source, 'chrompre', '')
        self.datatype = get_dict_value(source, 'datatype', '')
        self.buff = buff
        self.block_chrom = ""
        self.block_maxpos = 0
        self.block = {}

    def get_varkey_from_vcfrecord(self, vcfrecord):
        chrom = vcfrecord[0].replace('chr', '').strip()
        pos = int(vcfrecord[1])
        ref = vcfrecord[3].strip()
        alt = vcfrecord[4].strip()
        varkey = vcfrecord[1] + '_' + ref + '>' + alt
        return chrom, pos, ref, alt, varkey

    def get_annotmap(self, vcfrecord, varkeyset):
        chrom, pos, ref, alt, varkey = self.get_varkey_from_vcfrecord(vcfrecord)

        if chrom != self.block_chrom or pos > self.block_maxpos:
            print("load block", str(pos) + '~' + str(pos + MAXBUFF), self.sourcename, self.block_maxpos)
            self.load_annotblock(chrom, pos, varkeyset)

        try:
            annotmap = self.block[varkey]
        except KeyError:
            annotmap = {}
        return annotmap

    def load_annotblock(self, chrom, pos, varkeyset={}):
        self.block = {}
        if self.fileformat == "tab" or self.fileformat == "tsv" or self.fileformat == "bed":
            self.load_annotblock_tabformat(chrom, pos, varkeyset)
            # annotmap = self.load_annotdata_tabformat(chrom, pos, ref, alt)
        elif self.fileformat == "info":
            self.load_annotblock_infoformat(chrom, pos, varkeyset)
            # annotmap = self.load_annotdata_infoformat(chrom, pos, ref, alt)
        elif self.fileformat == "infovar":
            self.load_annotblock_infovarformat(chrom, pos, varkeyset)
            # annotmap = self.load_annotdata_infovarformat(chrom, pos, ref, alt)

    def load_annotblock_tabformat(self, chrom, pos, varkeyset={}):
        if self.fileformat == 'tsv':
            colidx = TSVCOLIDX
            sidx = 4
        elif self.fileformat == 'tab':
            colidx = VCFCOLIDX
            sidx = 5
        elif self.fileformat == 'bed':
            sidx = 0

        # recs = self.tabixpointer.query(self.chrompre+chrom, pos, pos+self.buff)
        recs = self.tabixpointer.query(self.chrompre + chrom, pos, pos + MAXBUFF)
        for r1 in recs:
            d = {}
            for k in range(len(self.header)):
                d[self.header[k]] = conv_infovalue(r1[k + sidx].strip())
            varkey = r1[colidx['POS']] + '_' + r1[colidx['REF']] + '>' + r1[colidx['ALT']]

            if len(varkeyset.keys()) > 0:
                try:
                    varkeyset[varkey]
                    try:
                        self.block[varkey].append(d)
                    except KeyError:
                        self.block[varkey] = [d]
                except KeyError:
                    pass
            else:
                try:
                    self.block[varkey].append(d)
                except KeyError:
                    self.block[varkey] = [d]

            # print(list(varkeyset.keys())[:3])
            # print(pos, varkey, len(self.block.keys()))
            if len(self.block.keys()) % 10 == 0:
                print(len(self.block.keys()), varkey)
            if len(self.block.keys()) >= self.buff:
                self.block_maxpos = int(r1[colidx['POS']])
                self.block_chrom = chrom
                break

    def load_annotblock_infovarformat(self, chrom, pos, varkeyset={}):
        try:
            # recs = self.tabixpointer.query(self.chrompre+chrom, pos, pos+self.buff)
            recs = self.tabixpointer.query(self.chrompre + chrom, pos, pos + MAXBUFF)
            for r1 in recs:
                d = {}
                for infofield in r1[VCFCOLIDX['INFO']].split(';'):
                    if '=' in infofield:
                        infoarr = infofield.split('=')
                        fieldid = infoarr[0]
                        fieldvalue = infoarr[1]
                        d[fieldid] = conv_infovalue(fieldvalue.strip())

                varkey = r1[VCFCOLIDX['POS']] + '_' + r1[VCFCOLIDX['REF']] + '>' + r1[VCFCOLIDX['ALT']]
                if len(varkeyset.keys()) > 0:
                    try:
                        varkeyset[varkey]
                        try:
                            self.block[varkey].append(d)
                        except KeyError:
                            self.block[varkey] = [d]
                    except KeyError:
                        pass
                else:
                    try:
                        self.block[varkey].append(d)
                    except KeyError:
                        self.block[varkey] = [d]

                if len(self.block.keys()) % 10 == 0:
                    print(len(self.block.keys()), varkey)

                if len(self.block.keys()) >= self.buff:
                    self.block_maxpos = int(r1[VCFCOLIDX['POS']])
                    self.block_chrom = chrom
                    break
        except tabix.TabixError:
            pass

    def load_annotdata_tabformat(self, chrom, pos, ref, alt):
        recs = self.tabixpointer.query(self.chrompre + chrom, pos, pos)
        sdata = []
        for r1 in recs:
            if self.fileformat == 'tsv':
                colidx = TSVCOLIDX
                sidx = 4
            elif self.fileformat == 'tab':
                colidx = VCFCOLIDX
                sidx = 5
            elif self.fileformat == 'bed':
                sidx = 0

            if self.datatype == "region" or r1[colidx['REF']] == ref and r1[colidx['ALT']] == alt:
                d = {}
                for k in range(len(self.header)):
                    d[self.header[k]] = conv_infovalue(r1[k + sidx].strip())
                sdata.append(d)
        return sdata

    def load_annotdata_infoformat(self, chrom, pos, ref, alt):
        recs = self.tabixpointer.query(self.chrompre + chrom, pos, pos)
        sdata = []

        target_fieldid = self.sourcename
        if '.' in self.sourcename:
            target_fieldid = self.sourcename.split('.')[1].strip()

        for r1 in recs:
            if r1[VCFCOLIDX['REF']] == self.ref and r1[VCFCOLIDX['ALT']] == self.alt:
                for infofield in r1[VCFCOLIDX['INFO']].split(';'):
                    infoarr = infofield.split('=')
                    fieldid = infoarr[0]
                    if target_fieldid == fieldid:
                        for subfield in infoarr[1].strip().split(','):
                            arr = subfield.split('|')
                            if self.sourcename == "SNPEFF.LOF" or self.sourcename == "SNPEFF.NMD":
                                arr[0] = arr[0][1:]         # remove '('
                                arr[-1] = arr[-1][:-1]      # remove ')'
                            d = {}
                            for k in range(len(self.header)):
                                d[self.header[k]] = conv_infovalue(arr[k].strip())
                            sdata.append(d)
        return sdata

    def load_annotdata_infovarformat(self, chrom, pos, ref, alt):
        sdata = []
        try:
            recs = self.tabixpointer.query(self.chrompre + chrom, pos, pos)
            # recs = self.tabixpointer.query(self.chrompre+chrom, pos, pos+10*1000*1000)
            for r1 in recs:
                if r1[VCFCOLIDX['REF']] == ref and r1[VCFCOLIDX['ALT']] == alt:
                    d = {}
                    for infofield in r1[VCFCOLIDX['INFO']].split(';'):
                        if '=' in infofield:
                            infoarr = infofield.split('=')
                            fieldid = infoarr[0]
                            fieldvalue = infoarr[1]
                            d[fieldid] = conv_infovalue(fieldvalue.strip())
                    sdata.append(d)
        except tabix.TabixError:
            pass
        return sdata


class AnnotMapBlock():
    subembedded_keymap = {}

    def __init__(self, datastruct, datafileinfo, buff):
        self.datafileinfo = datafileinfo
        self.datastruct = datastruct
        self.subembedded_keymap = {}
        self.set_annotblock(buff)

    def set_annotblock(self, buff):
        self.annotblock = {}
        ds = self.datastruct
        for s1 in ds['source']:
            tp = self.datafileinfo['tps'][s1['name']]
            header = self.datafileinfo['headers'][s1['name']]
            self.annotblock[s1['name']] = AnnotBlock(s1, tp, header, buff)

    def get_annotmap(self, vcfrecord, varkeyset={}):
        annotmap = {}
        ds = self.datastruct
        for s1 in ds['source']:
            annotmap[s1['name']] = self.annotblock[s1['name']].get_annotmap(vcfrecord, varkeyset)
        return annotmap

    def set_variantsample_annot(self):
        pass

    def set_transcript_annot(self):
        pass

    def set_gene_annot(self):
        pass

    def get_annotvcf_block(self, vcfblock):
        vcfblockcont = ""
        varkeyset = {}
        for i in range(len(vcfblock)):
            v1 = vcfblock[i]
            varkeyset[v1[1] + '_' + v1[3] + '>' + v1[4]] = 1

        for i in range(len(vcfblock)):
            vcfrecord = self.get_annotvcf_record(vcfblock[i], varkeyset)
            vcfblockcont += vcfrecord + '\n'
        return vcfblockcont

    def get_annotvcf_record(self, vcfrecord, varkeyset={}):
        vcfrecord = self.cleanup_infofield(vcfrecord)
        annotmap = self.get_annotmap(vcfrecord, varkeyset)
        if 'merged_one_field' in self.datastruct.keys() and self.datastruct['merged_one_field'] != '':
            vcfrecord = self.update_info_with_annotmap_in_one_field(vcfrecord, annotmap)
        else:
            vcfrecord = self.update_info_with_annotmap(vcfrecord, annotmap)
        return '\t'.join(vcfrecord)

    def cleanup_infofield(self, vcfrecord):
        if vcfrecord[VCFCOLIDX['INFO']] != '' and vcfrecord[VCFCOLIDX['INFO']] != '.':
            vcfrecord[VCFCOLIDX['INFO']] += ';'
        if vcfrecord[VCFCOLIDX['INFO']] == '.':
            vcfrecord[VCFCOLIDX['INFO']] = ''
        return vcfrecord

    def update_info_with_annotmap_in_one_field(self, vcfrecord, annotmap):
        merged_annot = {}
        multi_fields_sourcename = ""
        no_multi_fields_source = 0
        for s1 in self.datastruct['source']:
            fieldselection = get_dict_value(s1, 'fieldselection', '')
            annotlist = []
            if s1['name'] in annotmap.keys():
                for subfield in annotmap[s1['name']]:
                    subannotlist = []
                    if fieldselection == "all":
                        for f1 in subfield.keys():
                            subannotlist.append(subfield[f1].strip())
                    else:
                        for finfo in s1['fields']:
                            f1 = finfo['name']
                            try:
                                a1 = subfield[f1].strip()
                            except KeyError:
                                a1 = ''
                            subannotlist.append(a1)
                    annotlist.append('|'.join(subannotlist))
            merged_annot[s1['name']] = annotlist
            if multi_fields_sourcename == "":
                multi_fields_sourcename = s1['name']
            if len(annotlist) > 1:
                no_multi_fields_source += 1
                multi_fields_sourcename = s1['name']

        if no_multi_fields_source <= 1:
            merged_annotlist = []
            for annotfield in merged_annot[multi_fields_sourcename]:
                for s1name in merged_annot.keys():
                    if s1name != multi_fields_sourcename:
                        if len(merged_annot[s1name]) > 0:
                            annotfield += '|' + merged_annot[s1name][0]
                        else:
                            annotfield += '|'
                merged_annotlist.append(annotfield)
            vcfrecord[VCFCOLIDX['INFO']] += self.datastruct['merged_one_field'] + '=' + ','.join(merged_annotlist) + ';'
        return vcfrecord

    def update_info_with_annotmap(self, vcfrecord, annotmap):
        for s1 in self.datastruct['source']:
            flag_post = False
            fieldselection = get_dict_value(s1, 'fieldselection', '')
            subembedded = get_dict_value(s1, 'subembedded', '')
            annotlist = []
            annotlistmap = {}

            flag_make_subembedded_keymap = False
            if subembedded != '':
                try:
                    self.subembedded_keymap[subembedded]
                except KeyError:
                    self.subembedded_keymap[subembedded] = []
                    flag_make_subembedded_keymap = True

            if s1['name'] == "VEP":
                annotlist = self.get_vep_annotlist(fieldselection, s1, annotmap)
                flag_post = True
                # print("subembedded_keymap:", self.subembedded_keymap)
            elif s1['name'] in annotmap.keys():
                annotlist = []
                for subfield in annotmap[s1['name']]:
                    subannotlist = []
                    if fieldselection == "all":
                        for f1 in subfield.keys():
                            subannotlist.append(subfield[f1].strip())
                    else:
                        for finfo in s1['fields']:
                            if get_dict_value(finfo, 'is_available', True):
                                f1 = finfo['name']
                                try:
                                    a1 = subfield[f1].strip()
                                except KeyError:
                                    a1 = ''

                                if get_dict_value(finfo, 'is_list', False):
                                    a1 = "~".join(a1.split(conv_infovalue(finfo['delimiter'])))
                                subannotlist.append(a1)

                                if get_dict_value(finfo, 'subembedded_key', False):
                                    if flag_make_subembedded_keymap:
                                        self.subembedded_keymap[subembedded].append(a1)
                                    else:
                                        subembedded_key = a1
                                        print('self.subembedded_keymap_index:',
                                              self.subembedded_keymap[subembedded].index(a1))

                                if a1 != '':
                                    flag_post = True
                    if len(subannotlist) > 0:
                        annotlist.append('|'.join(subannotlist))
                        if subembedded != '' and not flag_make_subembedded_keymap:
                            annotlistmap[self.subembedded_keymap[subembedded].index(
                                subembedded_key)] = '|'.join(subannotlist)
                            subannotlist_len = len(subannotlist)

            # TODO: ##########################################
            # print (">>>>>>>>>>>>>>",s1['name'],s1['fields'],annotlist)
            if flag_post:
                sourcename = get_sourcename(s1)
                if subembedded != '' and not flag_make_subembedded_keymap:
                    annotlist = self.get_annotlist_from_annotlistmap_in_subembedded(
                        annotlistmap, subannotlist_len, subembedded)
                vcfrecord[VCFCOLIDX['INFO']] += sourcename + '=' + ','.join(annotlist) + ';'
            return vcfrecord

    def get_annotlist_from_annotlistmap_in_subembedded(self, annotlistmap, subannotlist_len, subembedded):
        annotlist = []

        for k in range(len(self.subembedded_keymap[subembedded])):
            try:
                v1 = annotlistmap[k]
            except KeyError:
                v1 = ""
                for j in range(subannotlist_len):
                    v1 += '|'
            annotlist.append(v1)
        return annotlist

    def get_vep_annotlist(self, fieldselection, s1, annotmap):
        global entrezmap, refseqmap
        subembedded = get_dict_value(s1, 'subembedded', '')
        annotlist = []
        if s1['name'] in annotmap.keys():
            for subfield in annotmap[s1['name']]:
                subannotlist = []
                if fieldselection == "all":
                    for f1 in subfield.keys():
                        subannotlist.append(subfield[f1].strip())
                else:
                    if subfield['Feature'][:len('ENST')] == "ENST":
                        try:
                            entrez_xref = entrezmap[subfield['Gene']]
                        except KeyError:
                            entrez_xref = ''
                        try:
                            refseq_xref = refseqmap[subfield['Feature']]
                        except KeyError:
                            refseq_xref = ''

                        subfield['Feature_ncbi'] = refseq_xref
                        subfield['Gene_ncbi'] = entrez_xref

                        for finfo in s1['fields']:
                            if get_dict_value(finfo, 'is_available', True):
                                f1 = finfo['name']
                                try:
                                    a1 = subfield[f1].strip()
                                except KeyError:
                                    a1 = ''

                                if get_dict_value(finfo, 'subembedded_key', False):
                                    self.subembedded_keymap[subembedded].append(a1)

                                if get_dict_value(finfo, 'is_list', False):
                                    a1 = "~".join(a1.split(conv_infovalue(finfo['delimiter'])))
                                subannotlist.append(a1)

                if len(subannotlist) > 0:
                    annotlist.append('|'.join(subannotlist))
        return annotlist


class AnnotMap():
    def __init__(self, datastruct, datafileinfo):
        self.datastruct = datastruct
        self.tabixpointer = None
        self.gzpointer = None
        self.tchrom = ""
        self.flag_gzsource = True
        self.prev_record = []

    def load_tabixpointer(self):
        self.tabixpointer = tabix.open(self.datastruct['sourcefile'].replace("#CHROM#", self.tchrom))

    def load_gzpointer(self):
        if self.datastruct['sourcefile'].endswith('.gz'):
            self.flag_gzsource = True
            self.gzpointer = file_util.gzopen(self.datastruct['sourcefile'].replace("#CHROM#", self.tchrom))
        else:
            self.flag_gzsource = False
            self.gzpointer = open(self.datastruct['sourcefile'].replace("#CHROM#", self.tchrom), 'r')

    def get_annotrecord_from_gzpoint(self):
        rec = []
        if len(self.prev_record) > 2:
            rec = self.prev_record
        else:
            while True:
                line = self.gzpointer.readline()
                if self.flag_gzsource:
                    line = line.decode('UTF-8')
                if line[0] != '#':
                    rec = line.split('\t')
                    rec[-1] = rec[-1].strip()
                    break
                if line.strip() == '':
                    break
        return rec

    def write_annotvcf_with_vcfblock(self, vcfblock, fp):
        stime = time.time()
        for b1 in vcfblock:
            chrom = b1[0].replace('chr', '')
            pos = int(b1[1])
            ref = b1[3]
            alt = b1[4]
            if chrom != self.tchrom:
                self.tchrom = chrom
                self.load_gzpointer()

            if b1[7] == '.':
                b1[7] = ""

            while True:
                r1 = self.get_annotrecord_from_gzpoint()
                apos = int(r1[1])
                if apos == pos:
                    if r1[3] == ref and r1[4] == alt:
                        if b1[7].strip() != '':
                            b1[7] += ';'
                        b1[7] += r1[5]
                        self.prev_record = []
                        break
                elif apos > pos:
                    break
                else:
                    self.prev_record = []
            fp.write('\t'.join(b1) + '\n')
        etime = time.time()
        elapsed = etime - stime
        print('processed...', len(vcfblock), "elapsed:", elapsed)
        return len(vcfblock)

    # slow but useful for sparse VCF
    def write_annotvcf_with_vcfblock_tabix(self, vcfblock, fp, is_rm_unannotated = False):
        for b1 in vcfblock:
            chrom = b1[0].replace('chr', '')
            pos = int(b1[1])
            ref = b1[3]
            alt = b1[4]
            if chrom != self.tchrom:
                self.tchrom = chrom
                self.load_tabixpointer()

            flag_add = False
            recs = self.tabixpointer.query(chrom, pos, pos+1)
            for r1 in recs:
                # print(r1,b1)
                if int(r1[1]) == pos and r1[3] == ref and r1[4] == alt:
                    if b1[7] == '.':
                        b1[7] = ""
                    if b1[7].strip() != '':
                        b1[7] += ';'
                    b1[7] += r1[5]
                    flag_add = True
                    break
            if flag_add or not is_rm_unannotated:
                fp.write('\t'.join(b1) + '\n')
        
        return len(vcfblock)


class VCFBlockReader():
    def __init__(self, vcf, blocksize=10000):
        self.vcf = vcf
        self.blocksize = blocksize
        self.fp = file_util.gzopen(self.vcf)
        self.is_gz = self.vcf.endswith('.gz')
        self.eof = False
        self.total_variant = 0
        self.i_variant = 0

    def get_header(self, add_header):
        headercont = ""
        for line in file_util.gzopen(self.vcf):
            if self.is_gz:
                line = line.decode('UTF-8')
            if line[0] == "#":
                if line[:len('#CHROM')] == "#CHROM":
                    headercont += add_header
                if line[:len('##contig=')] != '##contig=':
                    headercont += line
            elif line.strip() != '':
                self.total_variant += 1
        return headercont

    def get_block(self):
        block = []
        while True:
            line = self.fp.readline()
            if self.is_gz:
                line = line.decode('UTF-8')

            if line.strip() != '':
                if line[0] != '#':
                    arr = line.split('\t')
                    arr[-1] = arr[-1].strip()
                    if "," in arr[4]:
                        for a1 in vcf_util.split_multiallelic_variants(arr):
                            block.append(a1)
                            self.i_variant += 1
                    else:
                        block.append(arr)
                        self.i_variant += 1
            else:
                self.eof = True
                break
            if len(block) >= self.blocksize:
                break

        return block


class AnnotVCF():
    vcfheader = None
    datafileinfo = {}
    datastruct = {}
    opt = {}

    def __init__(self, opt):
        self.vcfheader = None
        self.datafileinfo = {}
        self.datafileinfo['files'] = {}
        self.datafileinfo['tps'] = {}
        self.datafileinfo['headers'] = {}

        self.opt = opt
        self.datastruct = file_util.load_json(opt['ds'])
        self.blocksize = opt['blocksize']

        self.set_datafile()

    def get_datafile_header(self, datafile, sourcename, dataformat):
        global DATAHEADER
        header = []
        try:
            header = DATAHEADER[sourcename]
        except KeyError:
            for line in file_util.gzopen(datafile):
                line = line.decode('UTF-8')
                if dataformat == "tab":
                    if line[:len("#CHR")].lower() == "#chr":
                        arr = line[1:].strip().split('\t')
                        header = arr[5:]
                        break
                elif dataformat == "tsv":
                    if line[:len("#CHR")].lower() == "#chr":
                        arr = line[1:].strip().split('\t')
                        header = arr[4:]
                        break
                elif dataformat == "bed":
                    if line[:len("#CHR")].lower() == "#chr":
                        arr = line[1:].strip().split('\t')
                        header = arr
                        break
                elif dataformat == "info":
                    if line[:len("##INFO=")] == "##INFO=":
                        infoID = line[len("##INFO="):].split(',')[0].replace('<ID=', '')
                        target_fieldid = sourcename
                        if '.' in sourcename:
                            target_fieldid = sourcename.split('.')[1].strip()
                        if infoID == target_fieldid:
                            headercont = line.split(
                                'Description="')[-1].split(':')[-1].strip().replace('">', '').replace("'", "")
                            if ' | ' in headercont:
                                header = headercont.split(' | ')
                            else:
                                header = headercont.split('|')
                            break
                elif dataformat == "infovar":
                    if line[:len("##INFO=")] == "##INFO=":
                        infoID = line[len("##INFO="):].split(',')[0].replace('<ID=', '')
                        header.append(infoID)
                if line[0] != '#':
                    break
        DATAHEADER[sourcename] = header
        return header

    def set_datafile(self, chrom='1'):
        for s1 in self.datastruct['source']:
            if 'datafile' in s1.keys():
                datafile = s1['datafile'].replace('#CHROM#', chrom)
                if file_util.is_exist(datafile):
                    tp = tabix.open(datafile)
                    header = self.get_datafile_header(datafile, s1['name'], s1['format'])
                    self.datafileinfo['files'][s1['name']] = datafile
                    self.datafileinfo['tps'][s1['name']] = tp
                    self.datafileinfo['headers'][s1['name']] = header
                else:
                    # print('Error: File not exist.', datafile)
                    pass

    def get_annot_header(self):
        cont = ""
        if 'merged_one_field' in self.datastruct.keys() and self.datastruct['merged_one_field'] != '':
            fields = []
            for s1 in self.datastruct['source']:
                if get_dict_value(s1, 'is_available', True):
                    fieldselection = get_dict_value(s1, 'fieldselection', '')
                    if fieldselection == "all":
                        if 'use_sourcename_in_fieldname' in self.datastruct.keys():
                            if not self.datastruct['use_sourcename_in_fieldname']:
                                fields.extend(DATAHEADER[s1['name']])
                        else:
                            for f1 in DATAHEADER[s1['name']]:
                                fields.append(s1['name'] + "_" + DATAHEADER[s1['name']][f1])

                    else:
                        for f1 in s1['fields']:
                            if 'use_sourcename_in_fieldname' in self.datastruct.keys() and \
                                    not self.datastruct['use_sourcename_in_fieldname']:
                                if 'name2' in f1.keys():
                                    fields.append(f1['name2'])
                                else:
                                    fields.append(f1['name'])
                            else:
                                if 'name2' in f1.keys():
                                    fields.append(s1['name'] + "_" + f1['name2'])
                                else:
                                    fields.append(s1['name'] + "_" + f1['name'])
            cont += "##INFO=<ID=" + self.datastruct['merged_one_field'] + ",Number=.,Type=String,Description=\"" + \
                    self.datastruct['merged_one_field_desc'] + ". Format:'" + '|'.join(fields) + "' \">\n"
        else:
            for s1 in self.datastruct['source']:
                if get_dict_value(s1, 'is_available', True):
                    fieldselection = get_dict_value(s1, 'fieldselection', '')
                    fields = []
                    if fieldselection == "all":
                        fields = DATAHEADER[s1['name']]
                        fieldlist = []
                        for fieldname in fields:
                            fieldlist.append({'name': fieldname, 'desc': ''})
                        s1['fields'] = fieldlist
                    else:
                        for f1 in s1['fields']:
                            if get_dict_value(f1, 'is_available', True):
                                if get_dict_value(f1, 'name2', False):
                                    fields.append(f1['name2'])
                                else:
                                    fields.append(f1['name'])
                    sourcename = get_sourcename(s1)
                    cont += "##INFO=<ID=" + sourcename + ",Number=.,Type=String"

                    cont += ",Description=\"" + s1['desc'] + ". "
                    if get_dict_value(s1, 'subembedded', '') != '':
                        cont += "Subembedded:'" + s1['subembedded'] + "':"
                    cont += "Format:'" + '|'.join(fields) + "' \">\n"
        return cont

    def get_version_info(self):
        ds = self.datastruct

        cont = ""
        cont += '##MUTANNO=<ID=MUTANNO'
        if 'version' in ds.keys():
            cont += ',Version="' + ds['version'] + '"'
        if 'version_date' in ds.keys():
            cont += ',Date="' + ds['version_date'] + '"'
        if 'data_version' in ds.keys():
            cont += ',DataVersion="' + ds['data_version'] + '"'
        if  'data_version_date' in ds.keys():
            cont += ',DataDate="' + ds['data_version_date'] + '"'
        cont += '>\n'

        for s1 in ds['source']:
            sourcename = get_sourcename(s1)
            cont += '##MUTANNO=<ID=' + sourcename 
            if 'version' in s1.keys():
                cont += ',Version="' + s1['version'] + '"'
            if 'version_date' in s1.keys():
                cont += ',Date="' + s1['version_date'] + '"'
            cont += '">\n'
        return cont

    def run_by_loading_source(self):
        stime = time.time()
        annotmap = AnnotMapBlock(self.datastruct, self.datafileinfo, self.blocksize)
        etime = time.time()
        elapsed = etime - stime
        print(cnt, elapsed, 'sec')

        '''
            varno += 1
            # if self.vcfheader is None:
            #     self.vcfheader = vcf_util.VCFHEADER(headercont)
            #     f.write(self.vcfheader.header.strip() + '\n')
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()

            chrom = arr[0].replace('chr', '')
            if prev_chrom != chrom:
                self.set_datafile(chrom)

            if flag_block:
                vcfblock.append(arr)
                if len(vcfblock) >= self.opt['buff']:
                    f.write(am.get_annotvcf_block(vcfblock))
                    vcfblock = []
            else:
                f.write(am.get_annotvcf_record(arr) + '\n')
            prev_chrom = chrom
            # print(arr[:5])

            if varno > 100:
                break

        if len(vcfblock) > 0:
            f.write(am.get_annotvcf_block(vcfblock))
            vcfblock = []
        '''

    def run(self):
        f = open(self.opt['out'], 'w')
        stime0 = time.time()
        total_varno = 0
        am = AnnotMap(self.datastruct, self.datafileinfo)
        vblock = VCFBlockReader(self.opt['vcf'], self.opt['blocksize'])
        f.write(vblock.get_header(self.get_version_info() + self.get_annot_header()))
        while(not vblock.eof):
            stime1 = time.time()
            # varno += am.write_annotvcf_with_vcfblock(vblock.get_block(), f)
            varno = am.write_annotvcf_with_vcfblock_tabix(vblock.get_block(), f, self.opt['remove_unannotated_variant'])
            total_varno += varno
            etime1 = time.time()
            elapsed1 = etime1 - stime1
            log = 'processed... ' + str(total_varno) + '/' + str(vblock.total_variant) 
            log += " elapsed:" + str(round(elapsed1,2)) + "s"
            # log += " time:" + str( (elapsed1/varno) * (vblock.total_variant - total_varno))
            print(log)

        etime0 = time.time()
        elapsed0 = etime0 - stime0
        print("total:", elapsed0, ", time per variant:", elapsed0 / total_varno,
              ", time for 4M:", str(round(elapsed0 / total_varno * 4000000 / 3600, 4)) + "hr")
        f.close()
        print("Saved " + self.opt['out'])


if __name__ == "__main__":
    pass
