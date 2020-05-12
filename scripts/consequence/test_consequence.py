#!/usr/bin/env python
# -*- coding: utf-8 -*-
# check_consequence.py
# made by Daniel Minseok Kwon
# 2020-04-28 15:04:23
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
import tabix
import mutanno

class ENSEMBL_GENEANNOTATION_GTF:
    # type_list = ['gene', 'transcript', 'exon', 'CDS', 'start_codon', 'stop_codon', 'five_prime_utr', 'three_prime_utr', 'Selenocysteine']
    def __init__(self, gtf):
        self.tp = tabix.open(gtf)

    
    def parse_fields(self, arr):
        d = {}
        d['chrom'] = arr[0].strip()
        d['spos'] = int(arr[3].strip())
        d['epos'] = int(arr[4].strip())
        d['type'] = arr[2].strip()
        d['strand'] = arr[6].strip()
        for f1 in arr[-1].split(';'):
            if f1.strip() != "":
                arr = f1.strip().split(' "')
                d[arr[0]] = arr[1].replace('"','')
        return d

    def get_data(self, chrom, spos, epos, margin=5000, include_all_types=True):
        datlist = {}
        
        for rec in self.tp.query(chrom, spos-margin, epos+margin):
            d = self.parse_fields(rec)
            try:
                datlist[d['type']]
            except KeyError:
                datlist[d['type']] = []
            datlist[d['type']].append(d)

        for ti in range(len(datlist['transcript'])):
            transcript = datlist['transcript'][ti]
            exonlist = []
            if transcript['spos'] < (spos - margin) or transcript['epos'] > (epos + margin):
                for rec in self.tp.query(chrom, transcript['spos'], transcript['epos']):
                    d = self.parse_fields(rec)
                    if d['type'] == 'exon' and d['transcript_id'] == transcript['transcript_id']:
                        exonlist.append(d)
            else:
                for exon in datlist['exon']:
                    if exon['transcript_id'] == datlist['transcript'][ti]['transcript_id']:
                        exonlist.append(exon)
            datlist['transcript'][ti]['exon'] = exonlist


            intronlist = []
            for ei in range(len(exonlist)-1):
                exon1 = exonlist[ei]
                exon2 = exonlist[ei+1]
                intron = {}
                intron['chrom'] = exon1['chrom']
                intron['spos'] = exon1['epos'] + 1
                intron['epos'] = exon2['spos'] - 1
                intron['strand'] = exon1['strand']
                intron['type'] = 'intron'
                intron['transcript_id'] = exon1['transcript_id']
                intronlist.append(intron)

            datlist['transcript'][ti]['intron'] = intronlist

        return datlist

class ENSEMBL_REGULATION_GFF:
    def __init__(self, gff):
        self.data = {}
        self.load_gff(gff)

    def load_gff(self, gff):
        i = 0
        for line in file_util.gzopen(gff):
            line = file_util.decodeb(line)
            if line[0] != '#':
                i += 1

                d = self.parse_fields(line)
                try:
                    self.data[d['chrom']]
                except KeyError:
                    self.data[d['chrom']] = []

                self.data[d['chrom']].append(d)

        # print(self.data.keys())
        # print('loading complete.')

    def parse_fields(self, line):
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()

        d = {}
        d['chrom'] = arr[0].strip()
        d['spos'] = int(arr[3].strip())
        d['epos'] = int(arr[4].strip())
        d['type'] = arr[2].strip()
        # d['score'] = float(arr[5].strip())
        d['strand'] = arr[6].strip()
        for f1 in arr[-1].split(';'):
            if f1.strip() != "":
                arr = f1.strip().split('=')
                d[arr[0]] = arr[1]
        d['ID'] = d['ID'].split(':')[1]
        return d

    def get_data(self, chrom, spos, epos, margin=0):
        margin = 0
        datlist = []
        for d in self.data[chrom]:
            if epos >= (d['spos']-margin) and spos <= (d['epos'] + margin):
                datlist.append(d)
            # if spos > (d['epos'] + margin):
            #     break
            
        return datlist

class ENSEMBL_MOTIF_GFF:
    def __init__(self, gff):
        self.data = {}
        self.tp = tabix.open(gff)

    def parse_fields(self, arr):
        d = {}
        d['chrom'] = arr[0].strip()
        d['spos'] = int(arr[3].strip())
        d['epos'] = int(arr[4].strip())
        d['type'] = arr[2].strip()
        d['score'] = float(arr[5].strip())
        d['strand'] = arr[6].strip()
        for f1 in arr[-1].split(';'):
            if f1.strip() != "":
                arr = f1.strip().split('=')
                d[arr[0]] = arr[1]
        return d

    def get_data(self, chrom, spos, epos):
        datlist = []
        records = self.tp.query(chrom, spos, epos)
        for rec in records:
            d = self.parse_fields(rec)
            datlist.append(d)
        return datlist




class MUTANNO_VARIANT_CONSEQUENCE:
    def __init__(self, ensgene_gtf_file, ens_regulation_gff="", ens_motif_gff=""):
        self.geneannot = ENSEMBL_GENEANNOTATION_GTF(ensgene_gtf_file)
        
        if ens_regulation_gff != "":
            self.has_regulation = True
            self.regulation = ENSEMBL_REGULATION_GFF(ens_regulation_gff)
        else:
            self.has_regulation = False

        if ens_motif_gff != "":
            self.has_motif = True
            self.motif = ENSEMBL_MOTIF_GFF(ens_motif_gff)
        else:
            self.has_motif = False
        

    def get_spos_epos(self, pos, ref, alt):
        spos = pos - 1
        if ref == '' or alt == '':
            epos = pos
        elif ref[0] == alt[0]:
            epos = spos + len(ref)
            spos += 1
        else:
            epos = spos + len(ref)
        return spos, epos

    def init_annotdata(self):
        annotdata = {'Gene':'', 'Feature':'', 'Feature_type':'', 'SYMBOL':'', 'Consequence':''}
        return annotdata

    def is_exonal_region(self, spos, epos, exon):
        flag = False
        # print(spos, epos, exon['spos'], exon['epos'])
        if spos < exon['epos'] and epos >= exon['spos']:
            flag = True
        return flag

    # A sequence variant in which a change has occurred within the region of the splice site, 
    # either within 1-3 bases of the exon or 3-8 bases of the intron
    def is_splice_region(self, spos, epos, exon, exon_number, total_exon_number):
        flag = False
        spos1 = spos + 1
        # print(int(exon['exon_number']))
        # print(spos1, epos, exon['spos'],exon['epos'], exon_number, total_exon_number)
        # print('exon:',exon)
        en = exon_number
        es = exon['strand']

        if (en != 1 and es == '+') or (en != total_exon_number and es == '-'):
            if spos1 >= (exon['spos'] - 8) and spos1 <= (exon['spos'] - 3):
                flag = True
            if spos1 >= (exon['spos']) and spos1 < (exon['spos'] + 3):
                flag = True

        if (en != total_exon_number and es == '+') or (en != 1 and es == '-'):
            if epos <= (exon['epos'] + 8) and epos >= (exon['epos'] + 3):
                flag = True
            
            if epos <= (exon['epos']) and epos > (exon['epos'] - 3):
                flag = True
        return flag


    # A splice variant that changes the 2 base region at the 3' end of an intron
    def is_splice_acceptor(self, spos, epos, exon, exon_number, total_exon_number):
        flag = False
        spos1 = spos + 1
        # print(spos, epos, exon_number,total_exon_number, exon)
        if exon['strand'] == '+':
            if exon_number != 1 and epos >= (exon['spos'] - 2) and spos1 <= (exon['spos'] - 1):
                flag = True
        elif exon['strand'] == '-':
            if exon_number != 1 and epos >= (exon['epos'] + 1) and spos1 <= (exon['epos'] + 2):
                flag = True
        return flag

    # A splice variant that changes the 2 base region at the 5' end of an intron
    def is_splice_donor(self, spos, epos, exon, exon_number, total_exon_number):
        flag = False
        spos1 = spos + 1
        # print("%%>", exon['strand'], exon_number, total_exon_number, spos, epos, exon['spos'], exon['epos'])
        if exon['strand'] == '+':
            if exon_number != total_exon_number and epos >= (exon['epos'] + 1) and spos1 <= (exon['epos'] + 2):
                flag = True
        elif exon['strand'] == '-':
            if exon_number != total_exon_number and epos >= (exon['spos'] - 2) and spos1 <= (exon['spos'] - 1):
                flag = True
        return flag

    def get_exon_consequence_tag(self, spos, epos, exonlist, gene_biotype, ctag):
        total_exon_number=len(exonlist)
        is_exon_variant = False
        for exon_i in range(total_exon_number):
            exon = exonlist[exon_i]
            # print(exon)

            if exon['strand'] == '+':
                exon_number = exon_i + 1
            elif exon['strand'] == '-':
                exon_number = total_exon_number - exon_i

            if self.is_exonal_region(spos, epos, exon):
                is_exon_variant = True

            if self.is_splice_region(spos, epos, exon, exon_number, total_exon_number):
                
                if 'splice_region_variant' not in ctag:
                    ctag.append('splice_region_variant')
            
            if self.is_splice_acceptor(spos, epos, exon, exon_number, total_exon_number):
                if 'splice_acceptor_variant' not in ctag:
                    ctag.append('splice_acceptor_variant')

            if self.is_splice_donor(spos, epos, exon, exon_number, total_exon_number):
                if 'splice_donor_variant' not in ctag:
                    ctag.append('splice_donor_variant')
        
        if is_exon_variant:
            if gene_biotype == 'protein_coding':
                ctag.append('exon_variant')
            else:
                ctag.append('non_coding_transcript_exon_variant')

        return ctag

    def is_intron_region(self, spos, epos, intron, include_slice = True):
        flag = False
        if include_slice and spos < intron['epos'] and epos >= intron['spos']:
            flag = True
        if not include_slice and spos < (intron['epos']-2) and epos >= (intron['spos']+2):
            flag = True
        return flag

    def get_intron_consequence_tag(self, spos, epos, intronlist, gene_biotype, ctag):
        total_intron_number=len(intronlist)
        is_intron_variant = False
        is_intron_variant_include_splice = False
        for intron_i in range(total_intron_number):
            intron = intronlist[intron_i]

            print('intron:', intron)
            if self.is_intron_region(spos, epos, intron, False):
                is_intron_variant = True
            if self.is_intron_region(spos, epos, intron, True):
                is_intron_variant_include_splice = True

        if is_intron_variant_include_splice:
            if is_intron_variant:
                ctag.append('intron_variant')
            if gene_biotype != 'protein_coding':
                ctag.append('non_coding_transcript_variant')

        return ctag


    def cal_consequence(self, spos, epos, vartype, transcript, margin):
        ctag = []
        spos1 = spos + 1

        print('########', spos, epos, transcript['spos'], transcript['epos'])

        if spos1 <= transcript['spos'] and epos >= transcript['epos']:
            if vartype == "DEL":
                ctag.append('transcript_ablation')
        elif spos1 <= transcript['epos'] and epos >= transcript['spos']:
            print('>>>>1')
            
            ctag = self.get_exon_consequence_tag(spos, epos, transcript['exon'], transcript['gene_biotype'], ctag)
            ctag = self.get_intron_consequence_tag(spos, epos, transcript['intron'], transcript['gene_biotype'], ctag)


        elif epos <= transcript['spos'] and spos1 >= (transcript['spos'] - margin):
            print('>>>>2')
            if transcript['strand'] == '+':
                ctag.append('upstream_gene_variant')
            elif transcript['strand'] == '-':
                ctag.append('downstream_gene_variant')
        elif epos >= transcript['epos'] and spos1 <= (transcript['epos'] + margin):
            print('>>>>3')
            if transcript['strand'] == '+':
                ctag.append('downstream_gene_variant')
            elif transcript['strand'] == '-':
                ctag.append('upstream_gene_variant')
        else:
            print('>>>>4')
            pass

        return '&'.join(ctag)

    def get_geneannot(self, chrom, spos, epos, vartype):
        margin = 5000
        annotlist = []
        geneannot_list = self.geneannot.get_data(chrom, spos, epos, margin, include_all_types=True)

        if 'transcript' in geneannot_list.keys():
            for transcript in geneannot_list['transcript']:

                annotdata = self.init_annotdata()
                annotdata['Gene'] = transcript['gene_id']
                annotdata['Feature'] = transcript['transcript_id']
                annotdata['Feature_type'] = 'Transcript'
                annotdata['SYMBOL'] = transcript['gene_name']
                annotdata['Biotype'] = transcript['transcript_biotype']
                annotdata['Gene_biotype'] = transcript['gene_biotype']

                annotdata['Consequence'] = self.cal_consequence(spos, epos, vartype, transcript, margin)

                annotdata['SPOS'] = transcript['spos']
                annotdata['EPOS'] = transcript['epos']
                if annotdata['Consequence'] != "":
                    annotlist.append(annotdata)
        else:
            annotdata = self.init_annotdata()
            annotdata['Consequence'] = "intergenic_variant"
            annotlist.append(annotdata)

        # print(self.ensgenelist[chrom])
        return annotlist

    def get_regulation(self, chrom, spos, epos):
        annotlist = []
        regulation_list = self.regulation.get_data(chrom, spos, epos)

        for regulation in regulation_list:
            annotdata = self.init_annotdata()
            annotdata['Feature'] = regulation['ID']
            annotdata['Feature_type'] = 'RegulatoryFeature'
            annotdata['Consequence'] = 'regulatory_region_variant'
            annotlist.append(annotdata)
        return annotlist

    def get_motif(self, chrom, spos, epos):
        annotlist = []
        motif_list = self.motif.get_data(chrom, spos, epos)

        for motif in motif_list:
            # print(motif)
            annotdata = self.init_annotdata()
            annotdata['Feature'] = motif['stable_id']
            annotdata['Feature_type'] = 'MotifFeature'
            annotdata['Consequence'] = motif['type'] + '_variant'
            annotlist.append(annotdata)
        return annotlist


    def get_vartype(self, ref, alt):
        lr = len(ref)
        la = len(alt)
        vartype = ""
        if lr == la and lr == 1:
            vartype = "SNV"
        elif lr > la:
            vartype = "DEL"
        elif lr < la:
            vartype = "INS"
        return vartype

    def get_consequence(self, chrom, pos, ref='', alt=''):
        spos, epos = self.get_spos_epos(pos, ref, alt)
        vartype = self.get_vartype(ref, alt)
        annotlist = self.get_geneannot(chrom, spos, epos, vartype)
        if self.has_regulation:
            annotlist.extend(self.get_regulation(chrom, spos, epos))
        if self.has_motif:
            annotlist.extend(self.get_motif(chrom, spos, epos))
        return annotlist


def is_diff_consequence(consequence1, consequence2):
    is_diff = False
    clist1 = consequence1.split('&')
    clist2 = consequence2.split('&')
    for c1 in clist1:
        if c1 not in clist2:
            is_diff = True
            break
    if not is_diff:
        for c2 in clist2:
            if c2 not in clist1:
                is_diff = True
                break
    return is_diff
        

def test_consequence(ensgene_gtf, ens_regulatory_gff, ens_motif_gff, vcf):
    # mt_consequence = MUTANNO_VARIANT_CONSEQUENCE(ensgene_gtf, ens_regulatory_gff, ens_motif_gff)
    mt_consequence = MUTANNO_VARIANT_CONSEQUENCE(ensgene_gtf)

    mt_seq_reader = mutanno.datasource.MUTANNO_DATA_SEQUENTIAL_READER(vcf)
    i = 0
    prev_pos = 0
    k = 0
    for mtdata in mt_seq_reader:
        i += 1
        # if prev_pos != mtdata['POS']:
        # if i % 100000 == 0:
        # if mtdata['POS'] > 12185:
        #     break
        # if mtdata['POS'] >= 12057 and mtdata['POS'] <= 12060:
        # if mtdata['POS'] == 12058:
        if True:
        # 946483
        # 1108052
        # 1205366
            flag = True
            ctag = mt_consequence.get_consequence(mtdata['CHROM'], mtdata['POS'], mtdata['REF'], mtdata['ALT'])

            # print(mtdata['VEP'])
            print("[",mtdata['CHROM'], mtdata['POS'], mtdata['REF'], mtdata['ALT'],"]")
            vep_transcript = {}
            vep_map = {}
            if 'VEP' in mtdata.keys():
                for v1 in mtdata['VEP']:
                    if v1['Feature_type'] == "Transcript" or v1['Feature_type'] == "":
                        vep_transcript[v1['Feature']] = v1['Consequence']
                        vep_map[v1['Feature']] = v1            

                mt_transcript = {}
                for ct1 in ctag:
                    if ct1['Feature_type'] == "Transcript" or v1['Feature_type'] == "":
                        mt_transcript[ct1['Feature']] = ct1['Consequence']
                        if ct1['Feature'] in vep_transcript.keys() and vep_transcript[ct1['Feature']] not in ["mature_miRNA_variant"]:
                            if is_diff_consequence(ct1['Consequence'], vep_transcript[ct1['Feature']]):
                                print('==> X', ct1['Feature'], vep_transcript[ct1['Feature']], ct1['Consequence'])
                                print(vep_map[ct1['Feature']])
                                print(ct1)
                                flag = False

                
                if len(vep_transcript.keys()) == len(mt_transcript.keys()) and flag:
                    # print('OK')
                    # print(vep_map)
                    # print(vep_transcript)
                    # print(mt_transcript)
                    print("\t".join(list(mt_transcript.values())))
                    pass
                else:
                    print ("=============")
                    print('ERROR')
                    print(len(vep_transcript.keys()), len(mt_transcript.keys()))
                    print ("#VEP#")
                    print(vep_transcript)
                    print ("-------------")
                    print ("#MUTANNO#")
                    print(mt_transcript)
                    print ("-------------")
                    print ()
                    break
                

            prev_pos = mtdata['POS']
            k += 1

            # break
            
        if k > 100000:
            break
        

if __name__ == "__main__":
    import proc_util
    import file_util
    # ensgene_gtf = "/home/mk446/mutanno/DATASOURCE/ENSEMBL/hg38/Homo_sapiens.GRCh38.99.gtf.gz"
    ensgene_gtf = "/home/mk446/mutanno/DATASOURCE/ENSEMBL/hg38/Homo_sapiens.GRCh38.99.sorted.gtf.gz"
    ens_regulation_gff = "/home/mk446/mutanno/DATASOURCE/ENSEMBL/hg38/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz"
    ens_motif_gff = "/home/mk446/mutanno/DATASOURCE/ENSEMBL/hg38/Homo_sapiens.GRCh38.motif_features.gff.gz"
    vcf = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/v0.4_200408_chrom/microannot_datasource.#CHROM#.v0.4_200408.tsi.gz"
    vcf = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/v0.4_200408/test.sorted.tsi.gz"
    test_consequence(ensgene_gtf, ens_regulation_gff, ens_motif_gff, vcf.replace('#CHROM#', '1'))
