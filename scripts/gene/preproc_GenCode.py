#!/usr/bin/env python
# -*- coding: utf-8 -*-
# preproc_GenCode.py
# made by Daniel Minseok Kwon
# 2020-03-31 10:16:35
#########################
import sys
import os
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path = "/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path = "/ms1/bin/python_lib"
else:
    sys_path = "/home/mk446/bin/python_lib"
sys.path.append(sys_path)
sys.path.append("..")


TAGLIST = "exon|CDS|start_codon|stop_codon|five_prime_UTR|three_prime_UTR".split('|')
# DELIMITER = "~"
DELIMITER = ";"

def info_pars(fieldcont):
    m = {}
    for field in fieldcont.split(';'):
        if '=' in field:
            arr2 = field.split('=')
            m[arr2[0].strip()] = arr2[1].strip()
    arr_id = m['ID'].split('.')
    m['id'] = arr_id[0].strip()
    m['id_version'] = arr_id[1].strip()
    return m

class ENST():
    def __init__(self, arr):
        info = info_pars(arr[8])
        arr_id = info['ID'].split('.')
        self.id = info['id']
        self.id_version = info['id_version']
        self.regions = {}
        self.pos = arr[3] + '-' + arr[4]
        self.refseq_list = []
        self.refseq_protein_list = []
        self.transcript_seq = ""
        self.protein_seq = ""
        self.ensp_idv = ""
 
    def add_line_annotation(self, arr):
        info = info_pars(arr[8])
        ftype = arr[2]
        pos = arr[3] + '-' + arr[4]
        try:
            self.regions[ftype].append(pos)
        except KeyError:
            self.regions[ftype] = [pos]
    def add_line_refseq(self, arr):
        self.refseq_list.append(arr[1])
        if len(arr) > 2 and arr[2].strip() != "":
            self.refseq_protein_list.append(arr[2])

    def get_region(self, ftype):
        global DELIMITER
        cont = []
        try:
            cont = self.regions[ftype]
        except KeyError:
            pass
        return DELIMITER.join(cont)

    def get_refseq(self):
        return DELIMITER.join(self.refseq_list)

    def get_refseq_protein(self):
        return DELIMITER.join(self.refseq_protein_list)


class ENSG():
    def __init__(self, arr):
        info = info_pars(arr[8])
        arr_id = info['ID'].split('.')
        # print(info)
        self.id = info['id']
        self.idv = info['ID']
        self.id_version = info['id_version']
        self.gene_name = info['gene_name']
        self.gene_type = info['gene_type']
        self.chrom = arr[0].replace('chr','')
        self.spos = arr[3]
        self.epos = arr[4]
        self.strand = arr[6]
        try:
            self.hgnc_id = info['hgnc_id'].replace('HGNC:','').strip()
        except KeyError:
            self.hgnc_id = ''
        self.transcripts = {}
    
    def add_line_annotation(self, arr):
        info = info_pars(arr[8])
        ftype = arr[2]
        pos = arr[3] + '-' + arr[4]

        if ftype == "transcript":
            self.transcripts[info['ID']] = ENST(arr)
        else:
            et = self.transcripts[info['Parent']]
            et.add_line_annotation(arr)

    def get_transcriptlist(self):
        global TAGLIST
        rst = []
        transcriptlist = list(self.transcripts.keys())
        for enstidv in transcriptlist:
            et = self.transcripts[enstidv]
            cont = [enstidv]
            cont.append(et.pos)
            for tag in TAGLIST:
                cont.append(et.get_region(tag))

            cont.append(et.get_refseq())
            cont.append(et.get_refseq_protein())
            cont.append(et.ensp_idv)
            if et.transcript_seq == '':
                cont.append('')
            else:
                cont.append(str(len(et.transcript_seq)))
            
            cont.append(et.transcript_seq)
            if et.transcript_seq == '':
                cont.append('')
            else:
                cont.append(str(len(et.protein_seq)))
            cont.append(et.protein_seq)

            rst.append(cont)
        return rst


class GenCode():
    def __init__(self, path, version):
        self.path = path
        self.version = version
        self.ensglist = {}
        self.enstlist = {}
        
    
    def add_line_annotation(self, arr):
        ftype = arr[2]
        info = info_pars(arr[8])
        if ftype == "gene":
            if info['gene_id'] not in self.ensglist.keys():
                self.ensglist[info['gene_id']] = ENSG(arr)
        else:
            if ftype not in ["stop_codon_redefined_as_selenocysteine"]:
                eg = self.ensglist[info['gene_id']]
                eg.add_line_annotation(arr)
            # print(ftype)
    def mapping_transcript(self):
        for ensg_idv in self.ensglist.keys():
            eg = self.ensglist[ensg_idv]
            for enst_idv in eg.transcripts.keys():
                self.enstlist[enst_idv] = eg.transcripts[enst_idv]

    def add_line_refseq(self, arr):
        print(arr)
        enst_idv = arr[0]
        try:
            et = self.enstlist[enst_idv]
            et.add_line_refseq(arr)
        except KeyError:
            print('KeyError:', enst_idv)
            pass
    
    def save_header(self, f, output_type='gene'):
        global TAGLIST
        cont = []
        cont.append("#CHROM")
        cont.append("SPOS")
        cont.append("EPOS")
        cont.append("STRAND")
        cont.append("ensgid")
        cont.append("ensgid_version")
        cont.append("gene_name")
        cont.append("gene_type")
        cont.append("hgnc_id")

        theader = ["enstid","pos"]
        theader.extend(TAGLIST)
        theader.append("refseq")
        theader.append("refseq_protein")
        theader.append("enspid")
        theader.append("transcript_length")
        theader.append("transcript_sequence")
        theader.append("protein_length")
        theader.append("protein_sequence")

        if output_type == 'gene':
            cont.append('|'.join(theader))
        else:
            cont.extend(theader)
        f.write('\t'.join(cont) + '\n')

    def save_tsv(self, f):
        self.save_header(f, 'gene')
        for engs_idv in self.ensglist.keys():
            eg = self.ensglist[engs_idv]
            cont = []
            cont.append(eg.chrom)
            cont.append(eg.spos)
            cont.append(eg.epos)
            cont.append(eg.strand)
            cont.append(eg.id)
            cont.append(eg.idv)
            cont.append(eg.gene_name)
            cont.append(eg.gene_type)
            cont.append(eg.hgnc_id)

            cont_transcript = []
            transcriptlist = eg.get_transcriptlist()
            for transcript in transcriptlist:
                cont_transcript.append('|'.join(transcript))
            cont.append(','.join(cont_transcript))
            # print(eg.regions)
            f.write('\t'.join(cont) + '\n')

    def save_tsv_transcript_based(self, f):
        self.save_header(f, 'transcript')
        for engs_idv in self.ensglist.keys():
            eg = self.ensglist[engs_idv]
            transcriptlist = eg.get_transcriptlist()
            for transcript in transcriptlist:
                cont = []
                cont.append(eg.chrom)
                cont.append(eg.spos)
                cont.append(eg.epos)
                cont.append(eg.strand)
                cont.append(eg.id)
                cont.append(eg.idv)
                cont.append(eg.gene_name)
                cont.append(eg.gene_type)
                cont.append(eg.hgnc_id)
                cont.extend(transcript)
                f.write('\t'.join(cont) + '\n')

    def add_transcript_seq(self, enst_idv, seq):
        try:
            et = self.enstlist[enst_idv]
            et.transcript_seq = seq
        except KeyError:
            pass

    def add_protein_seq(self, enst_idv,ensp_idv, seq):
        try:
            et = self.enstlist[enst_idv]
            et.ensp_idv = ensp_idv
            et.protein_seq = seq
        except KeyError:
            pass


    def add(self, fname, addtype):
        if addtype == "transcript_seq":
            seq = ""
            enst_idv= ""
            for line in file_util.gzopen(fname):
                line = line.decode('UTF-8').strip()
                if line[0] == '>':
                    if len(seq) > 0 and enst_idv  != "":
                        self.add_transcript_seq(enst_idv, seq)
                        seq = ""
                        enst_idv = ""

                    enst_idv = line[1:].split('|')[0]
                else:
                    seq += line

            if len(seq) > 0 and enst_idv  != "":
                self.add_transcript_seq(enst_idv, seq)
        elif addtype == "protein_seq":
            seq = ""
            enst_idv= ""
            for line in file_util.gzopen(fname):
                line = line.decode('UTF-8').strip()
                if line[0] == '>':
                    if len(seq) > 0 and enst_idv  != "":
                        self.add_protein_seq(enst_idv,ensp_idv, seq)
                        seq = ""
                        enst_idv = ""
                    arr = line[1:].split('|')
                    ensp_idv = arr[0]
                    enst_idv = arr[1]
                else:
                    seq += line

            if len(seq) > 0 and enst_idv  != "":
                self.add_protein_seq(enst_idv,ensp_idv, seq)

        else:
            i = 0
            for line in file_util.gzopen(fname):
                line = line.decode('UTF-8')
                if line[0] != '#':
                    i += 1
                    arr = line.split('\t')
                    arr[-1] = arr[-1].strip()
                    if addtype == "annotation":
                        self.add_line_annotation(arr)
                    if addtype == "refseq":
                        self.add_line_refseq(arr)
                    
                    if i > 200:
                        pass
                        # break

def check_output(out):
    for line in open(out):
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        if line[0] == "#":
            infoheader = arr[-1].split('|')
        else:
            for f1 in arr[-1].split(','):
                info = f1.split('|')
                if len(info) != len(infoheader):
                    print(len(info), len(infoheader))
                    print(arr)
                    break

def preproc_GenCode(version, out, out_transcript_based):

    gc = GenCode(path, version)

    f = open(out, 'w')
    ft = open(out_transcript_based, 'w')

    gc.add(path + 'gencode.'+version+'.annotation.gff3.gz', 'annotation')
    gc.mapping_transcript()
    gc.add(path + 'gencode.'+version+'.metadata.RefSeq.gz', 'refseq')
    gc.add(path + 'gencode.'+version+'.transcripts.fa.gz', 'transcript_seq')
    gc.add(path + 'gencode.'+version+'.pc_translations.fa.gz', 'protein_seq')

    gc.save_tsv(f)
    gc.save_tsv_transcript_based(ft)

    f.close()
    ft.close()
    print('Saved', out)
    print('Saved', out_transcript_based)
    

if __name__ == "__main__":
    import preproc_util
    import proc_util
    import file_util
    path = preproc_util.DATASOURCEPATH + "/GENE/GENCODE/"
    version = "v33"
    version_date = "12_13_2019"
    out = path + "GenCode." + version + "."+version_date+".bed"
    out_transcript_based = path + "GenCodeTranscript." + version + "."+version_date+".bed"
    preproc_GenCode(version, out, out_transcript_based)
    # check_output(out, out_transcript_based)
