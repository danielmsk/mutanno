#!/usr/bin/env python
# -*- coding: utf-8 -*-
# mk_known_indel_vcf.py
# made by Min-Seok Kwon
# 2020-01-13 19:32:02
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


def mk_known_indel_vcf(f, vcf, chromidx, posidx, refidx, altidx):
    print("processing...", vcf)
    for line in file_util.gzopen(vcf):
        line = line.decode('UTF-8')
        if line[0] != '#':
            arr = line.split('\t')
            if not (len(arr[refidx]) == 1 and len(arr[altidx]) == 1):
                cont = 'chr' + arr[chromidx].replace('chr', '') + '\t' + arr[posidx] + \
                    '\t.\t' + arr[refidx] + '\t' + arr[altidx]
                cont += addfield
                f.write(cont + '\n')
            if ',' in arr[refidx] or ',' in arr[altidx]:
                print('Error:', arr, vcf)


if __name__ == "__main__":
    import file_util
    import seq_util
    out = "/home/mk446/mutanno/DATASOURCE/KNOWN_INDEL/hg38/known_indel2.vcf"

    addfield = '\t.\t.\t.\tGT:AD:DP\t0/1:0,30:60'
    f = open(out, 'w')
    header = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1'
    # f.write(header + '\n')

    vcf = "/home/mk446/mutanno/DATASOURCE/PATHOGENICITY/CADD/hg38/v1.5/ALL.TOPMed_freeze5_hg38_dbSNP.tsv.gz"
    mk_known_indel_vcf(f, vcf, 0, 1, 2, 3)

    vcf = "/home/mk446/mutanno/DATASOURCE/PATHOGENICITY/CADD/hg38/v1.5/InDels.tsv.gz"
    mk_known_indel_vcf(f, vcf, 0, 1, 2, 3)

    vcf = "/home/mk446/mutanno/DATASOURCE/POPULATION_AF/gnomAD/hg38/gnomad.genomes.r3.0.sites.AF_popmax.tsv.gz"
    mk_known_indel_vcf(f, vcf, 0, 1, 2, 3)

    vcf = "/home/mk446/mutanno/DATASOURCE/POPULATION_AF/gnomAD/hg38/r2.1.1/"
    vcf += "gnomad.exomes.r2.1.1.sites.liftover_grch38.AF_popmax.tsv.gz"
    mk_known_indel_vcf(f, vcf, 0, 1, 2, 3)

    vcf = "/home/mk446/mutanno/DATASOURCE/POPULATION_AF/UK10K/hg38/UK10K_COHORT.20160215.sites.liftover.sorted.vcf.gz"
    mk_known_indel_vcf(f, vcf, 0, 1, 3, 4)

    vcf = "/home/mk446/mutanno/DATASOURCE/POPULATION_AF/TOPMED/hg38/bravo-dbsnp-all.vcf.gz"
    mk_known_indel_vcf(f, vcf, 0, 1, 3, 4)

    # TODO: Split multiple variants in ESP6500
    # vcf = "/home/mk446/mutanno/DATASOURCE/POPULATION_AF/ESP6500/hg38/"
    # vcf += "ESP6500SI-V2-SSA137.GRCh38-liftover.hg38pos.sorted.vcf.gz"
    # mk_known_indel_vcf(f, vcf, 0, 1, 3, 4)

    # TODO: Split multiple variants in 1000G
    # vcf = "/home/mk446/mutanno/DATASOURCE/POPULATION_AF/1000GP/hg38/ALL.chr#CHROM#_GRCh38_sites.20170504.vcf.gz"
    # for chrom in seq_util.MAIN_CHROM_LIST:
    #     if chrom != 'MT':
    #         mk_known_indel_vcf(f, vcf.replace('#CHROM#', chrom), 0, 1, 3, 4)

    vcf = "/home/mk446/mutanno/DATASOURCE/POPULATION_AF/KAVIAR/hg38/"
    vcf += "Kaviar-160204-Public-hg38-trim.splited.sorted.vcf.gz"
    mk_known_indel_vcf(f, vcf, 0, 1, 3, 4)

    f.close()
