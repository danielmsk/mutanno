from pyliftover import LiftOver

CHROM_LIST = {}
CHROM_LEN = {}
CHROM_TYPE = {}
CHROM_TYPE_LIST = ['autosome', 'sex', 'mt', 'unlocalized', 'unplaced', 'haplotype', 'hla', 'ebv', 'decoy']
REF_SEQ_VERSION_LIST = ['b37', 'hg19', 'b37d5', 'b38', 'hg38', 'b38d']
MAIN_CHROM_LIST = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "M"]

CHROM_LEN['hg38'] = {}
CHROM_LEN['hg38']['1'] = 248956422
CHROM_LEN['hg38']['2'] = 242193529
CHROM_LEN['hg38']['3'] = 198295559
CHROM_LEN['hg38']['4'] = 190214555
CHROM_LEN['hg38']['5'] = 181538259
CHROM_LEN['hg38']['6'] = 170805979
CHROM_LEN['hg38']['7'] = 159345973
CHROM_LEN['hg38']['8'] = 145138636
CHROM_LEN['hg38']['9'] = 138394717
CHROM_LEN['hg38']['10'] = 133797422
CHROM_LEN['hg38']['11'] = 135086622
CHROM_LEN['hg38']['12'] = 133275309
CHROM_LEN['hg38']['13'] = 114364328
CHROM_LEN['hg38']['14'] = 107043718
CHROM_LEN['hg38']['15'] = 101991189
CHROM_LEN['hg38']['16'] = 90338345
CHROM_LEN['hg38']['17'] = 83257441
CHROM_LEN['hg38']['18'] = 80373285
CHROM_LEN['hg38']['19'] = 58617616
CHROM_LEN['hg38']['20'] = 64444167
CHROM_LEN['hg38']['21'] = 46709983
CHROM_LEN['hg38']['22'] = 50818468
CHROM_LEN['hg38']['X'] = 156040895
CHROM_LEN['hg38']['Y'] = 57227415


codon = {}
codon['UUU'] = "F"
codon['UUC'] = "F"
codon['UUA'] = "L"
codon['UUG'] = "L"
codon['CUU'] = "L"
codon['CUC'] = "L"
codon['CUA'] = "L"
codon['CUG'] = "L"
codon['AUU'] = "I"
codon['AUC'] = "I"
codon['AUA'] = "I"
codon['AUG'] = "M"
codon['GUU'] = "V"
codon['GUC'] = "V"
codon['GUA'] = "V"
codon['GUG'] = "V"
codon['UCU'] = "S"
codon['UCC'] = "S"
codon['UCA'] = "S"
codon['UCG'] = "S"
codon['CCU'] = "P"
codon['CCC'] = "P"
codon['CCA'] = "P"
codon['CCG'] = "P"
codon['ACU'] = "T"
codon['ACC'] = "T"
codon['ACA'] = "T"
codon['ACG'] = "T"
codon['GCU'] = "A"
codon['GCC'] = "A"
codon['GCA'] = "A"
codon['GCG'] = "A"
codon['UAU'] = "Y"
codon['UAC'] = "Y"
codon['UAA'] = "*"
codon['UAG'] = "*"
codon['CAU'] = "H"
codon['CAC'] = "H"
codon['CAA'] = "Q"
codon['CAG'] = "Q"
codon['AAU'] = "N"
codon['AAC'] = "N"
codon['AAA'] = "K"
codon['AAG'] = "K"
codon['GAU'] = "D"
codon['GAC'] = "D"
codon['GAA'] = "E"
codon['GAG'] = "E"
codon['UGU'] = "C"
codon['UGC'] = "C"
codon['UGA'] = "*"
codon['UGG'] = "W"
codon['CGU'] = "R"
codon['CGC'] = "R"
codon['CGA'] = "R"
codon['CGG'] = "R"
codon['AGU'] = "S"
codon['AGC'] = "S"
codon['AGA'] = "R"
codon['AGG'] = "R"
codon['GGU'] = "G"
codon['GGC'] = "G"
codon['GGA'] = "G"
codon['GGG'] = "G"

complementary = {}
complementary['A'] = 'T'
complementary['T'] = 'A'
complementary['G'] = 'C'
complementary['C'] = 'G'

BASE = ["A", "T", "G", "C"]
COMPBASE = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

# SIGNATURE
SINGLEBASE_PATTERN6 = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
SINGLEBASE_PATTERN = ['C>A', 'G>T', 'C>G', 'G>C', 'C>T', 'G>A', 'T>A', 'A>T', 'T>C', 'A>G', 'T>G', 'A>C']
SINGLEBASE_COMPLEMENTARY = {'G>T': 'C>A', 'C>A': 'C>A', 'C>G': 'C>G', 'G>C': 'C>G', 'C>T': 'C>T', 'G>A': 'C>T',
                            'T>A': 'T>A', 'A>T': 'T>A', 'T>C': 'T>C', 'A>G': 'T>C', 'T>C': 'T>C', 'A>G': 'T>C', 'T>G': 'T>G', 'A>C': 'T>G'}

EXOM_SEQ_COVERED_REGION_NimbleGen_SeqCap = 64190747
EXOM_SEQ_COVERED_REGION = EXOM_SEQ_COVERED_REGION_NimbleGen_SeqCap

# REF_FASTA = {}
# REF_FASTA['b37d5'] = "/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta"
# REF_FASTA['hg38'] = "/home/mk446/BiO/Install/GATK-bundle/2.8/hg38/Homo_sapiens_assembly38.fasta"
# REF_FASTA['hg19'] = "/home/mk446/BiO/Install/GATK-bundle/2.8/hg19/ucsc.hg19.fasta"
# REF_FASTA['b38d'] = "/home/mk446/BiO/Install/GATK-bundle/2.8/b38/GRCh38_full_analysis_set_plus_decoy_hla.fa"


def load_liftover(chainfile = ""):
    if chainfile == "":
        lo = LiftOver('hg38', 'hg19')
    else:
        lo = LiftOver(chainfile)
    return lo

def convert_coordinate(lo_hg38_hg19, chrom, pos):
    if 'chr' not in chrom:
        chrom = 'chr' + chrom
    if type(pos) != int:
        pos = int(pos)
    return lo_hg38_hg19.convert_coordinate(chrom, pos)



def load_refseq_info(seq_ver=''):
    global CHROM_LIST, CHROM_LEN, CHROM_TYPE, REF_SEQ_VERSION_LIST
    path = '/home/mk446/bin/python_lib/'
    if seq_ver == '':
        for seq_ver in REF_SEQ_VERSION_LIST:
            load_refseq_info(seq_ver)
    else:
        CHROM_LIST[seq_ver] = []
        CHROM_LEN[seq_ver] = {}
        CHROM_TYPE[seq_ver] = {}
        for line in open(path+'ref_'+seq_ver+'.txt'):
            arr = line.strip().split('\t')
            if arr[0] != "CHROM" and arr[0] != "":
                CHROM_LIST[seq_ver].append(arr[0])
                CHROM_LEN[seq_ver][arr[0]] = int(arr[1])
                CHROM_TYPE[seq_ver][arr[0]] = arr[2]

# load_refseq_info('b37d5')
# load_refseq_info('b38d')


def get_split_region(chunksize=1000000, seqver='b38d', tchrom= '', tspos=-9, tepos=-9):
    rst = []
    if tchrom != '' and tspos > 0:
        if tepos == -9:
            load_refseq_info(seqver)
            tepos = CHROM_LEN[seqver][tchrom]
        
        flag = True
        spos = tspos
        k2 = 0
        while flag:
            epos = spos + chunksize - 1
            if epos > tepos:
                epos = tepos
                flag = False
            k2 += 1
            rst.append((tchrom, spos, epos, k2))
            spos = epos + 1

    else:
        load_refseq_info(seqver)
        k2 = 1
        for chrom in CHROM_LIST[seqver]:
            if len(chrom) < 3:
                no_split = int(1.0 * CHROM_LEN[seqver][chrom] / chunksize)
                for k in range(1, no_split + 2):
                    spos = chunksize * (k - 1) + 1
                    epos = chunksize * (k)
                    if epos > CHROM_LEN[seqver][chrom]:
                        epos = CHROM_LEN[seqver][chrom]
                    rst.append((chrom, spos, epos, k2))
                    k2 += 1
    return rst


def codeAA(dnaseq3):
    global codon
    aa = ""
    if len(dnaseq3) == 3:
        aa = codon[dnaseq3.upper().replace('T', 'U')]
    return aa


def complementarySeq(seq, direction='+'):
    global complementary
    cseq = ""
    if direction == '-':
        for i in range(0, len(seq)):
            cseq += complementary[seq[i].upper()]
    else:
        for i in range(len(seq)-1, -1, -1):
            cseq += complementary[seq[i].upper()]

    return cseq

# DNA to Protein


def translate(dnaseq, orient="+"):
    global codon
    p = ""
    if orient == "-":
        r_dnaseq = dnaseq[::-1]

    for i in range(0, len(dnaseq), 3):
        if orient == "-":
            p += codeAA(complementarySeq(r_dnaseq[i:(i+3)]))
        else:
            p += codeAA(dnaseq[i:(i+3)])

    return p


# DNA to Protein
def translateAll(dnaseq):
    global codon
    p = {}
    p['0'] = ""
    p['1'] = ""
    p['2'] = ""
    p['-0'] = ""
    p['-1'] = ""
    p['-2'] = ""
    r_dnaseq = dnaseq[::-1]
    for i in range(0, len(dnaseq), 3):
        p['0'] += codeAA(dnaseq[i:(i+3)])
        p['1'] += codeAA(dnaseq[(i+1):(i+1+3)])
        p['2'] += codeAA(dnaseq[(i+2):(i+2+3)])
        p['-0'] += codeAA(complementarySeq(r_dnaseq[i:(i+3)]))
        p['-1'] += codeAA(complementarySeq(r_dnaseq[(i+1):(i+1+3)]))
        p['-2'] += codeAA(complementarySeq(r_dnaseq[(i+2):(i+2+3)]))

    return p


def get_chroSeqLen(seqlen_file="/ms2/refseq/hg19/chr/chr_len.txt"):
    m = {}
    for line in open(seqlen_file):
        arr = line.strip().split("\t")
        m[arr[0]] = int(arr[1])
    return m

# OUTDATED.... USE get_refseq()
# def get_seq(chrom, spos, epos, strand='+', seq_path="/home/mk446/BiO/Install/GATK-bundle/2.8/b37"):
#    fasta = seq_path + "/human_g1k_v37_decoy." + chrom + ".fasta"
#    #fasta = maestro.FASTA['b37_chrom'].replace("#CHROM#",chrom)
#    f = open(fasta, "r")
#    f.seek(spos-1) ### start 1
#    rst = f.read(epos-spos+1)
#    if strand == '-':
#        rst = complementarySeq(rst)
#    f.close()
#    return rst

# UCSC hg19 - outdated


def get_dnaseq(chro, spos, epos, orient="+", seq_path="/ms2/refseq/hg19/chr/", seqlen_file="/ms2/refseq/hg19/chr/chr_len.txt", report_pos=False):
    seq_file = seq_path + "chr" + chro + ".fa"

    if orient == "-":
        seqlen = get_chroSeqLen(seqlen_file)
        sspos = seqlen['chr'+chro] - epos
        eepos = seqlen['chr'+chro] - spos
    else:
        sspos = spos
        eepos = epos

    # print "sspos", sspos
    # print "eepos", eepos

    i = 0
    this_line_spos = 0
    sidx = 0
    eidx = 0
    seq = ""
    for line in open(seq_file):
        if i > 0:
            next_pos = this_line_spos + len(line.strip())
            if sspos >= this_line_spos and sspos < next_pos:
                sidx = sspos-this_line_spos
            if sidx > 0:
                seq += line.strip()
            if eepos >= this_line_spos and eepos < next_pos:
                eidx = eepos-this_line_spos
                break
            this_line_spos += len(line.strip())
        i += 1

    dnaseq = seq[sidx:(eepos-sspos)+sidx]
    # print seq, epos-sspos+1, len(dnaseq)
    # print sidx, eidx
    if (report_pos):
        return (dnaseq, sspos, eepos)
    else:
        return dnaseq


def get_base_info_from_bam(bam, chrom, pos1, ref="", alt=""):
    import pysam
    pysamAF = pysam.AlignmentFile(bam, "rb")
    return get_base_info_from_pysamAF(pysamAF, chrom, pos1, ref="", alt="")


def get_base_info_from_pysamAF(pysamAF, chrom, pos1, ref="", alt=""):
    pos = pos1-1
    #pos = pos1+1
    info = {}
    for x in pysamAF.pileup(chrom, pos-5, pos+5, min_base_quality=0):
        if x.reference_pos == pos:
            # print x.nsegments, x.reference_pos, len(x.pileups)
            info['total_depth'] = len(x.pileups)
            info['chrom'] = chrom
            info['pos'] = x.reference_pos
            info['base_info'] = []
            base_stat = get_base_info_from_bam_init_base_stat()
            base_stat_mapq20 = get_base_info_from_bam_init_base_stat()
            base_stat_baseq20 = get_base_info_from_bam_init_base_stat()
            base_stat_mapq20baseq20 = get_base_info_from_bam_init_base_stat()
            for pr in x.pileups:
                a = pr.alignment
                # print pr.query_position, a.query_sequence[pr.query_position], a.query_qualities[pr.query_position], a.mapq, a.is_reverse
                try:
                    base = {}
                    base['pos_in_read'] = pr.query_position
                    base['base'] = a.query_sequence[pr.query_position]
                    base['base_qual'] = a.query_qualities[pr.query_position]
                    base['MAPQ'] = a.mapq
                    base['is_reverse'] = a.is_reverse
                    info['base_info'].append(base)

                    base_stat[a.query_sequence[pr.query_position]] += 1
                    if base['MAPQ'] >= 20:
                        base_stat_mapq20[a.query_sequence[pr.query_position]] += 1
                    if base['base_qual'] >= 20:
                        base_stat_baseq20[a.query_sequence[pr.query_position]] += 1
                    if base['base_qual'] >= 20 and base['MAPQ'] >= 20:
                        base_stat_mapq20baseq20[a.query_sequence[pr.query_position]] += 1

                except TypeError:
                    pass
            info['base_stat'] = add_af_dp_from_base_stat(base_stat, ref, alt)
            info['base_stat_mapq20'] = add_af_dp_from_base_stat(base_stat_mapq20, ref, alt)
            info['base_stat_baseq20'] = add_af_dp_from_base_stat(base_stat_baseq20, ref, alt)
            info['base_stat_mapq20baseq20'] = add_af_dp_from_base_stat(base_stat_mapq20baseq20, ref, alt)
    return info


def get_baserange_info_from_bam(bam, chrom, spos, epos, flag_read_info=False):
    import pysam
    infolist = {}
    samfile = pysam.AlignmentFile(bam, "rb")
    ref = ''
    alt = ''
    for x in samfile.pileup(chrom, spos-50, epos+50, min_base_quality=0):
        if x.reference_pos >= spos-1 and x.reference_pos <= epos+1:
            # print x.nsegments, x.reference_pos, len(x.pileups)
            info = {}
            info['total_depth'] = len(x.pileups)
            info['chrom'] = chrom
            info['pos'] = x.reference_pos
            info['base_info'] = []
            base_stat = get_base_info_from_bam_init_base_stat()
            base_stat_mapq20 = get_base_info_from_bam_init_base_stat()
            base_stat_baseq20 = get_base_info_from_bam_init_base_stat()
            base_stat_mapq20baseq20 = get_base_info_from_bam_init_base_stat()
            for pr in x.pileups:
                a = pr.alignment
                # print pr.query_position, a.query_sequence[pr.query_position], a.query_qualities[pr.query_position], a.mapq, a.is_reverse
                try:
                    base = {}
                    base['pos_in_read'] = pr.query_position
                    base['base'] = a.query_sequence[pr.query_position]
                    base['base_qual'] = a.query_qualities[pr.query_position]
                    base['MAPQ'] = a.mapq
                    base['is_reverse'] = a.is_reverse
                    if flag_read_info:
                        info['base_info'].append(base)

                    base_stat[a.query_sequence[pr.query_position]] += 1
                    if base['MAPQ'] >= 20:
                        base_stat_mapq20[a.query_sequence[pr.query_position]] += 1
                    if base['base_qual'] >= 20:
                        base_stat_baseq20[a.query_sequence[pr.query_position]] += 1
                    if base['base_qual'] >= 20 and base['MAPQ'] >= 20:
                        base_stat_mapq20baseq20[a.query_sequence[pr.query_position]] += 1

                except TypeError:
                    pass
            info['base_stat'] = add_af_dp_from_base_stat(base_stat, ref, alt)
            info['base_stat_mapq20'] = add_af_dp_from_base_stat(base_stat_mapq20, ref, alt)
            info['base_stat_baseq20'] = add_af_dp_from_base_stat(base_stat_baseq20, ref, alt)
            info['base_stat_mapq20baseq20'] = add_af_dp_from_base_stat(base_stat_mapq20baseq20, ref, alt)
            infolist[x.reference_pos+1] = info
    return infolist


def add_af_dp_from_base_stat(base_stat, ref, alt):
    if ref != "":
        base_stat['DP'] = base_stat['A']+base_stat['T']+base_stat['G']+base_stat['C']+base_stat['N']
        if base_stat['DP'] > 0:
            base_stat['AF'] = base_stat[alt]/base_stat['DP']
        else:
            base_stat['AF'] = 0
        base_stat['AD'] = base_stat[alt]
    return base_stat


def get_base_info_from_bam_init_base_stat():
    base_stat = {}
    base_stat['A'] = 0
    base_stat['T'] = 0
    base_stat['G'] = 0
    base_stat['C'] = 0
    base_stat['N'] = 0
    return base_stat

## samf : samf = pysam.AlignmentFile(bam, "rb")
## chrom : arr[0]
## pos = int(arr[1])
## ref = arr[4]
## alt = arr[5]
# return giv : odds ratio


def get_giv(samf, chrom, pos, ref, alt):
    d = {}
    cnt_reads = load_reads(samf, chrom, pos)
    cnt_r1_ref = get_value(cnt_reads['r1'], ref)
    cnt_r1_alt = get_value(cnt_reads['r1'], alt)
    cnt_r2_ref = get_value(cnt_reads['r2'], ref)
    cnt_r2_alt = get_value(cnt_reads['r2'], alt)
    if cnt_r1_ref == 0:
        iv_r1 = 1.0
    else:
        iv_r1 = (1.0*cnt_r1_alt/cnt_r1_ref)
    if cnt_r2_ref == 0:
        iv_r2 = 1.0
    else:
        iv_r2 = (1.0*cnt_r2_alt/cnt_r2_ref)

    if iv_r2 == 0:
        giv = 1.0
    else:
        giv = iv_r1 / iv_r2

    d['cnt_reads'] = cnt_reads
    d['cnt_r1_ref'] = cnt_r1_ref
    d['cnt_r1_alt'] = cnt_r1_alt
    d['cnt_r2_ref'] = cnt_r2_ref
    d['cnt_r2_alt'] = cnt_r2_alt
    d['iv_r1'] = iv_r1
    d['iv_r2'] = iv_r2
    d['giv'] = giv
    return d


def get_value(m, k1):
    try:
        v1 = m[k1]
    except KeyError:
        v1 = 0
    return v1


def load_reads(samf, chrom, pos_plus_1):
    pos = pos_plus_1 - 1
    m = {}
    m['r1'] = {}
    m['r2'] = {}
    for x in samf.pileup(chrom, pos-1, pos+1):
        if x.reference_pos == pos:
            # print x.nsegments, x.reference_pos, len(x.pileups)
            for pr in x.pileups:
                a = pr.alignment
                if pr.query_position != None:
                    base = a.query_sequence[pr.query_position]
                    # print pr.query_position, a.query_sequence[pr.query_position], a.query_sequence, a.is_reverse, a.is_read1, a.is_read2
                    # print pr.query_position, a.query_sequence[pr.query_position], a.query_qualities[pr.query_position], a.mapq, a.is_reverse
                    if a.mapq >= 20 and not a.is_duplicate and not a.is_unmapped:
                        if a.is_read1:
                            try:
                                m['r1'][base] += 1
                            except KeyError:
                                m['r1'][base] = 1
                        elif a.is_read2:
                            try:
                                m['r2'][base] += 1
                            except KeyError:
                                m['r2'][base] = 1

    return m


# def get_refseq(chrom, spos, epos, seqver="b38d"):
#     from pyfasta import Fasta
#     #refseq = seq_path + "/human_g1k_v37_decoy.fasta"
#     spos = spos - 1
#     refseq = REF_FASTA[seqver]
#     f = Fasta(refseq)
#     #print (sorted(f.keys()))
#     chromlist = list(f.keys())
#     #print (chromlist[:30])
#     fasta_chrom = ""
#     for c1 in chromlist:
#         if c1.split(' ')[0] == chrom:
#             fasta_chrom = c1
#             break
#     #print (fasta_chrom)
#     #fasta_chrom = chrom + ' dna:chromosome chromosome:GRCh37:'+chrom+':1:'+str(CHROM_LEN['b37d5'][chrom])+':1'
#     refseq = f[fasta_chrom][spos:epos+1]
#     return refseq


############################################
# SEQ QC
############################################
def get_qcfile_info(fname, qctool):
    if qctool == 'stat':
        info = get_bamstat_info(fname)
    if qctool == 'cmm':
        info = get_cmm_info(fname)
    if qctool == 'hist':
        info = get_hist_info(fname)
    return info


def get_cmm_info(fname):
    # fname : .cmm.alignment_summary_metrics file
    info = {}

    for line in open(fname):
        if line[:len('FIRST_OF_PAIR')] == 'FIRST_OF_PAIR':
            arr = line.split('\t')
            info['first_of_pair'] = arr[1].strip()
        if line[:len('SECOND_OF_PAIR')] == 'SECOND_OF_PAIR':
            arr = line.split('\t')
            info['second_of_pair'] = arr[1].strip()
        if line[:len('PAIR')] == 'PAIR':
            arr = line.split('\t')
            info['total_reads'] = arr[1].strip()
            info['mean_read_length'] = arr[15].strip()

    info['genome_coverage'] = str(round(int(info['total_reads']) * int(info['mean_read_length']) / CHROM_TOTAL_LEN, 3))
    info['exom_coverage'] = str(round(int(info['total_reads']) * int(info['mean_read_length']) / EXOM_SEQ_COVERED_REGION, 3))

    insert_size_metrics_header = []
    flag_metrics = False
    flag_size_dist = False
    info['insert_size_dist'] = ""
    for line in open(fname.replace('.cmm.alignment_summary_metrics', '.cmm.insert_size_metrics')):
        if flag_metrics:
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            for k in range(len(arr)):
                info[insert_size_metrics_header[k]] = arr[k].strip()
            flag_metrics = False
            break
        if flag_size_dist:
            line = line.strip()
            if line(line) > 2:
                info['insert_size_dist'] += line+'\n'

        if line[:len('MEDIAN_INSERT_SIZE')] == 'MEDIAN_INSERT_SIZE':
            insert_size_metrics_header = line.split('\t')
            insert_size_metrics_header[-1] = insert_size_metrics_header[-1].strip()
            flag_metrics = True

        if line[:len('insert_size')] == 'insert_size':
            flag_size_dist = True

    ### TODO .cmm.quality_by_cycle_metrics
    ### TODO .cmm.quality_distribution_metrics
    ### TODO .cmm.base_distribution_by_cycle_metrics
    ### TODO .cmm.quality_distribution_metrics
    return info


def get_hist_info(fname):
    info = {}
    read_depth10 = 0
    read_depth30 = 0
    read_depth50 = 0
    read_depth70 = 0
    read_depth100 = 0
    read_depth200 = 0
    for line in open(fname):
        if line[:len('genome')] == 'genome':

            arr = line.strip().split('\t')
            read_depth = int(arr[1])
            covered_region = arr[2]
            # print (arr,read_depth,read_depth)
            if read_depth == 0:
                info['uncovered_region'] = covered_region
            elif read_depth <= 400:
                info['covered_region (dp='+str(read_depth)+')'] = covered_region
            if read_depth > 10:
                read_depth10 += int(covered_region)
            if read_depth > 30:
                read_depth30 += int(covered_region)
            if read_depth > 50:
                read_depth50 += int(covered_region)
            if read_depth > 70:
                read_depth70 += int(covered_region)
            if read_depth > 100:
                read_depth100 += int(covered_region)
            if read_depth > 200:
                read_depth200 += int(covered_region)
    info['covered_region (dp>10)'] = str(read_depth10)
    info['covered_region (dp>30)'] = str(read_depth30)
    info['covered_region (dp>50)'] = str(read_depth50)
    info['covered_region (dp>70)'] = str(read_depth70)
    info['covered_region (dp>100)'] = str(read_depth100)
    info['covered_region (dp>200)'] = str(read_depth200)

    return info


def get_bamstat_info(fname):
    info = {}
    flag_all_chrom = False
    flag_only_main = False
    flag_xy = False
    seqversion = "b37"

    for line in open(fname):
        if '===ALL CHROM==' in line:
            flag_all_chrom = True

        if '=ONLY MAIN CHROM (1~22,X,Y,MT)==' in line:
            flag_only_main = True
            flag_all_chrom = False

        if '==X,Y CHROM==' in line:
            flag_xy = True
            flag_only_main = False
            flag_all_chrom = False

        if "READ_LEN" in line:
            info['read_length'] = line.strip().replace('READ_LEN:', '').strip()

        if not flag_only_main and not flag_all_chrom and not flag_xy:
            arr = line.split('\t')
            # print (arr)
            if arr[0] != "CHROM":
                info[arr[0] + '_chrom_mapped_reads'] = arr[2].replace(',', '')
                info[arr[0] + '_chrom_unmapped_reads'] = arr[3].replace(',', '')
                info[arr[0] + '_chrom_total_reads'] = arr[4].replace(',', '')
                if arr[0] == "1" and arr[1] == "248,956,422":
                    seqversion = "b38"
                if arr[0] == "1" and arr[1] == "249,250,621":
                    seqversion = "b37"
                if arr[0] == "hs37d5":
                    seqversion = "b37d5"
                if arr[0] == "6_apd_hap1":
                    seqversion = "hg19"
                if arr[0] == "19_KI270938v1_alt":
                    seqversion = "hg38"
                if arr[0] == "Un_JTFH01001998v1_decoy":
                    seqversion = "b38d"

                # if line[0] == "X":
                #     info['X_chrom_total_reads'] = arr[4].replace(',','')
                # if line[0] == "Y":
                #     info['Y_chrom_total_reads'] = arr[4].replace(',','')

        if flag_all_chrom and "COVERAGE:" in line:
            info['all_chrom_coverage'] = line.strip().replace('COVERAGE:', '').strip()
        if flag_all_chrom and "TOTAL READS:" in line:
            info['all_chrom_total_reads'] = line.strip().replace('TOTAL READS:', '').replace(',', '').strip()
        if flag_all_chrom and "MAPPED:" in line and not "UNMAPPED:" in line:
            info['all_chrom_mapped_reads'] = line.strip().replace('MAPPED:', '').replace(',', '').strip()
        if flag_all_chrom and "UNMAPPED:" in line:
            info['all_chrom_unmapped_reads'] = line.strip().replace('UNMAPPED:', '').replace(',', '').strip()

        if flag_only_main and "COVERAGE:" in line:
            info['main_chrom_coverage'] = line.strip().replace('COVERAGE:', '').split('(')[0].strip()
        if flag_only_main and "TOTAL READS:" in line:
            info['main_chrom_total_reads'] = line.strip().replace('TOTAL READS:', '').replace(',', '').split('(')[0].strip()
        if flag_only_main and "MAPPED:" in line and not "UNMAPPED:" in line:
            info['main_chrom_mapped_reads'] = line.strip().replace('MAPPED:', '').replace(',', '').split('(')[0].strip()
        if flag_only_main and "UNMAPPED:" in line:
            info['main_chrom_unmapped_reads'] = line.strip().replace('UNMAPPED:', '').replace(',', '').split('(')[0].strip()

        if flag_xy and "X/Y RATIO:" in line:
            info['xy_ratio'] = line.strip().replace('X/Y RATIO:', '').replace(',', '').split('(')[0].strip()
        if flag_xy and "EXT. GENDER:" in line:
            info['ext_gender'] = line.strip().replace('EXT. GENDER:', '').replace(',', '').split('(')[0].strip()

    info['all_chrom_mapped_ratio'] = str(round(int(info['all_chrom_mapped_reads'])/int(info['all_chrom_total_reads'])*100, 4))+'%'
    if int(info['main_chrom_total_reads']) > 0:
        info['main_chrom_mapped_ratio'] = str(round(int(info['main_chrom_mapped_reads'])/int(info['main_chrom_total_reads'])*100, 4))+'%'
    else:
        info['main_chrom_mapped_ratio'] = "0"
    info['main_chrom_read_ratio'] = str(round(int(info['main_chrom_total_reads'])/int(info['all_chrom_total_reads'])*100, 4))+'%'
    info['seqversion'] = seqversion
    return info
