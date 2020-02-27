#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s04_check_missvar_in_tsv.py
# made by Min-Seok Kwon
# 2020-01-15 15:08:39
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


def get_split_region(fastaobj, chunksize=1000000):
    chunklist = []
    k2 = 1
    for chrom in list(fastaobj.keys()):
        if len(chrom) < 6:
            no_split = int(1.0 * len(fastaobj[chrom]) / chunksize)
            for k in range(1, no_split + 2):
                spos = chunksize * (k - 1) + 1
                epos = chunksize * (k)
                if epos > len(fastaobj[chrom]):
                    epos = len(fastaobj[chrom])
                chunklist.append((chrom, spos, epos, k2))
                k2 += 1
    return chunklist


def get_inputvcf_path(chrom, spos, epos):
    out = os.path.join(path, chrom, str(int(spos / 1000000)),
                       chrom + '_' + str(spos) + '_' + str(epos) + '.vcf')
    out = os.path.abspath(out)
    return out


def s04_check_missvar_in_tsv(chrom, spos, epos):
    fastaobj = Fasta(fastafile, as_raw=True, sequence_always_upper=True)
    inputvcf = get_inputvcf_path(chrom, spos, epos)
    tsvgz = inputvcf.replace('.vcf', '.tsv.gz')

    print(inputvcf)
    print('loading..', tsvgz)
    varmap = {}
    if file_util.is_exist(tsvgz):
        for line in file_util.gzopen(tsvgz):
            line = line.decode('UTF-8')
            if line[0] != '#':
                arr = line.split('\t')
                varmap[arr[1] + "_" + arr[2] + "_" + arr[3]] = 1
    print('#var', len(varmap.keys()))

    contall = header
    missflag = False
    for k in range(spos, epos + 1):
        ref = fastaobj[chrom][k - 1]
        if ref != 'N' and ref.strip() != '':
            cont = chrom
            cont += '\t' + str(k)
            cont += '\t' + '.'
            cont += '\t' + ref
            cont += '\t#ALT#'
            cont += '\t' + '.'
            cont += '\t' + '.'
            cont += '\t' + '.'
            cont += '\t' + 'GT:AD:DP'
            cont += '\t' + '0/1:0,30:60'
            cont += '\n'
            for alt in ['A', 'T', 'G', 'C']:
                if ref != alt:
                    contall += cont.replace('#ALT#', alt)
                    try:
                        tmp = varmap[str(k) + '_' + ref + '_' + alt]
                    except KeyError:
                        missflag = True
    if missflag:
        file_util.fileSave(inputvcf, contall, 'w')
        print('Saved..', inputvcf)
    else:
        print('No missing..')
        file_util.fileSave(tsvgz + '.nomiss', '', 'w')


def run():
    fastaobj = Fasta(fastafile, as_raw=True, sequence_always_upper=True)
    chunklist = get_split_region(fastaobj, chunksize)
    for (chrom, spos, epos, k2) in chunklist:
        cmd = "python /home/mk446/mutanno/SRC/scripts/precal_vep/s04_check_missvar_in_tsv.py " + \
            chrom + ' ' + str(spos) + ' ' + str(epos)
        print(cmd)


if __name__ == "__main__":
    import file_util
    chunksize = 100000
    header = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
    path = "/home/mk446/mutanno/PRECALVEP/"
    fastafile = "/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh38d1/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    if len(sys.argv) == 1:
        run()
    else:
        s04_check_missvar_in_tsv(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
