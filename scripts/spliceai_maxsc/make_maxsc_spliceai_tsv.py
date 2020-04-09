#!/usr/bin/env python
# -*- coding: utf-8 -*-
# make_maxsc_spliceai_tsv.py
# made by Min-Seok Kwon
# 2020-01-10 04:44:05
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


def make_maxsc_spliceai_tsv_from_rawvcf(raw_vcf):
    out = raw_vcf.replace('.vcf.gz', '') + '.maxsc.tsv'
    f = open(out, 'w')
    header = '#CHROM\tPOS\tREF\tALT\tMaxScore'
    f.write(header + '\n')
    i = 0
    for line in file_util.gzopen(raw_vcf):
        line = line.decode('UTF-8')
        if line[0] != '#':
            i += 1
            arr = line.split('\t')
            arr2 = arr[7].split('|')
            maxsc = ''
            maxscf = -1
            for sc in arr2[2:6]:
                if maxscf < float(sc):
                    maxsc = sc
            if maxsc != "0.00":
                cont = arr[0] + '\t' + arr[1] + '\t' + arr[3] + '\t' + arr[4] + '\t' + maxsc
                f.write(cont + '\n')
            # break
            if i % 1000000 == 0:
                print(i, cont)
    f.close()


def run_chrom():
    for chrom in seq_util.MAIN_CHROM_LIST:
        if chrom == "MT":
            chrom = "M"
        cmd = "python /home/mk446/mutanno/SRC/scripts/spliceai_maxsc/make_maxsc_spliceai_tsv.py " + chrom
        print(cmd)


def make_maxsc_spliceai_tsv_from_chromtsv(raw_file, chrom):
    out = raw_file.replace('#CHROM#', chrom).replace('.tsv.gz', '') + '.maxsc.tsv'
    raw_file = raw_file.replace('#CHROM#', chrom)
    f = open(out, 'w')
    i = 0
    header = []
    prev_cont = ""
    prev_pos = ""
    prev_maxscf = 0
    for line in file_util.gzopen(raw_file):
        line = line.decode('UTF-8')
        if line[0] == '#':
            header = line.split('\t')
            f.write(line.strip() + '\tMAXDS\tMAXDSNAME\n')
        else:
            i += 1
            arr = line.split('\t')
            pos1 = '_'.join(arr[:4])
            maxsc = ''
            maxscname = ''
            maxscf = -1
            k = 4
            for sc in arr[5:9]:
                k += 1
                scf = float(sc)
                if maxscf < scf:
                    maxsc = sc
                    maxscf = scf
                    maxscname = header[k]
                elif maxscf == scf:
                    maxscname += ',' + header[k]
            if maxsc == '0.00':
                maxscname = ''

            
            cont = line.strip() + '\t' + maxsc + '\t' + maxscname
            if pos1 == prev_pos:
                if prev_maxscf < maxscf:
                    prev_cont = cont
                    prev_maxscf = maxscf
                    prev_pos = pos1
            else:
                # cont = line.strip() + '\t' + maxsc + '\t' + maxscname
                if prev_cont != "":
                    f.write(prev_cont + '\n')
                    
                prev_cont = cont
                prev_maxscf = maxscf
                prev_pos = pos1
            
            # break
            # if i % 1000 == 0:
            if i % 1000000 == 0:
                print(i, cont)
                # break

    f.write(prev_cont + '\n')
    f.close()
    print('Saved', out)


if __name__ == "__main__":
    import file_util
    import seq_util
    # raw_vcf = "/home/mk446/mutanno/DATASOURCE/SPLICING/SpliceAI/hg38/spliceai_scores.raw.snv.hg38.vcf.gz"
    # make_maxsc_spliceai_tsv_from_rawvcf(raw_vcf)

    # run_chrom()
    raw_file = "/home/mk446/mutanno/DATASOURCE/SPLICING/SpliceAI/hg38/spliceai.20191004.hg38.#CHROM#.tsv.gz"
    make_maxsc_spliceai_tsv_from_chromtsv(raw_file, sys.argv[1])
