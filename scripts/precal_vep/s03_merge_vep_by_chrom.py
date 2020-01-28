#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s03_merge_vep_by_chrom.py
# made by Min-Seok Kwon
# 2020-01-14 18:18:44
#########################

import sys
import os
from pyfaidx import Fasta
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path = "/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path = "/ms1/bin/python_lib"
else:
    sys_path = "/home/mk446/bin/python_lib"
sys.path.append(sys_path)


def s03_merge_vep_by_chrom(chrom):
    runsh = 's03_chrom' + chrom + '.sh'
    print('Saving..', runsh)
    f = open(runsh, 'w')
    out = path + "vep.98.hg38." + chrom + ".tsv"
    i = 0
    for k in range(400):
        vcfmap = {}
        for tsv in file_util.walk(path + "chr" + chrom + "/" + str(k) + "/", '.tsv.vcf.gz'):
            k1 = int(tsv.split('/')[-1].split('_')[1])
            vcfmap[tsv] = k1
        (ks, vs) = struct_util.sortdict(vcfmap)

        for tsv in ks:
            if i == 0:
                cmd = "zcat " + tsv + " > " + out + ";"
            else:
                cmd = "zcat " + tsv + " | grep -v '^#' >> " + out + ";"
            i += 1
            f.write(cmd + '\n')
    cmd = "sleep 120;"
    f.write(cmd + '\n')
    cmd = 'tabixgz ' + out
    # cmd = "tabix -f -p vcf " + tsv
    f.write(cmd + '\n')
    f.close()


def s03_merge_vep_by_chrom_pipe_bgzip(chrom):
    logfile = path + "vep.98.hg38." + chrom + ".tsv.gz.log"
    file_util.fileSave(logfile, '', 'w')
    i = 0
    # fasta = Fasta(FASTA, as_raw=True, sequence_always_upper=True)
    # ref = fasta[chrom][k - 1]
    # if ref != 'N' and ref.strip() != '':

    for k in range(400):
        vcfmap = {}
        for tsv in file_util.walk(path + "chr" + chrom + "/" + str(k) + "/", '.vep.txt.gz'):
            k1 = int(tsv.split('/')[-1].split('_')[1])
            vcfmap[tsv] = k1
        (ks, vs) = struct_util.sortdict(vcfmap)

        for tsv in ks:
            j = 0
            log = str(i + 1) + ': ' + tsv
            file_util.fileSave(logfile, log + '\n', 'a')
            for line in file_util.gzopen(tsv):
                line = line.decode('UTF-8')
                if line[:len('#CHROM')] == "#CHROM" or line[0] != '#':
                    arr = line.split('\t')
                    arr[-1] = arr[-1].strip()
                    line = arr[0].replace('chr', '') + '\t' + '\t'.join(arr[1:5]) + '\t' + arr[7][4:] + '\n'
                if (i == 0) or (i > 0 and line[0] != '#'):
                    print(line, end='')
                j += 1
                if j % 100000 == 0:
                    file_util.fileSave(logfile, '\t' + line[:30] + '\n', 'a')
            i += 1

    file_util.fileSave(logfile, 'Done.\n', 'a')


def run():
    for chrom in seq_util.MAIN_CHROM_LIST:
        if chrom == "MT":
            chrom = "M"
        outgz = path + "vep.99.hg38." + chrom + ".tsi.gz"
        cmd = "python /home/mk446/mutanno/SRC/scripts/precal_vep/s03_merge_vep_by_chrom.py " + chrom
        cmd += " | bgzip -c > " + outgz + ";"
        cmd += "sleep 120;"
        cmd += "tabix -f -p vcf " + outgz + ";"
        runsh = 's03_chrom' + chrom + '.sh'
        file_util.fileSave(runsh, cmd + '\n', 'w')
        print('Saved..', runsh)


if __name__ == "__main__":
    import struct_util
    import file_util
    import seq_util

    path = "/home/mk446/mutanno/PRECALVEP/"
    spath = "/home/mk446/mutanno/SRC/scripts/precal_vep/"
    FASTA = "/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh38d1/GRCh38_full_analysis_set_plus_decoy_hla.fa"

    # for chrom in seq_util.MAIN_CHROM_LIST:
    #     if chrom == "MT":
    #         chrom = "M"
    #     s03_merge_vep_by_chrom(chrom)

    if len(sys.argv) == 1:
        run()
    else:
        s03_merge_vep_by_chrom_pipe_bgzip(sys.argv[1])
