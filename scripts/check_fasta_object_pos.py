#!/usr/bin/env python
# -*- coding: utf-8 -*-
# check_fasta_object_pos.py
# made by Min-Seok Kwon
# 2020-01-13 21:33:52
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


def check_fasta_object_pos(fasta, cadd):
    # print(fasta.keys())
    for line in file_util.gzopen(cadd):
        line = line.decode('UTF-8')
        if line[0] != '#':
            arr = line.split('\t')
            chrom = arr[0]
            pos = int(arr[1])
            ref = arr[2]
            fref = fasta['chr' + chrom][pos - 1]
            if ref != fref:
                print(pos, ref, fref)


if __name__ == "__main__":
    import proc_util
    import file_util
    fasta_file = "/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh38d1/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    fasta = Fasta(fasta_file, as_raw=True, sequence_always_upper=True)
    cadd = "/home/mk446/mutanno/DATASOURCE/PATHOGENICITY/CADD/hg38/v1.5/whole_genome_SNVs.tsv.gz"
    check_fasta_object_pos(fasta, cadd)
