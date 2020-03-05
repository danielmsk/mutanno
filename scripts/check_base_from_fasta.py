#!/usr/bin/env python
# -*- coding: utf-8 -*-
# check_base.py
# made by Daniel Minseok Kwon
# 2020-03-03 14:07:04
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
from pyfaidx import Fasta
import tabix


def check_base(tsvgz, baseidx):
    

    fobj = Fasta(fasta, as_raw=True, sequence_always_upper=True)

    # tb = tabix.open(tsvgz)
    # recs = tb.querys("Y:1-1250000")
    # recs = tb.querys("1:1-1250000")
    prev = ""
    f = open(tsvgz + '.basecheck','w')
    # for r1 in recs:
    i = 0
    for line in file_util.gzopen(tsvgz):
        line = line.decode('UTF-8')
        if line[0] != '#':
            r1 = line.split('\t')
            r1[-1] = r1[-1].strip()
            if r1[1] != prev:
                i += 1
                chrom = 'chr'+r1[0]
                # chrom = r1[0]
                ref = fobj[chrom][int(r1[1])-1]
                if r1[baseidx] != ref:
                    cont = [r1[0], r1[1], r1[baseidx], ref]
                    f.write('\t'.join(cont) + '\n')
                    # break
                if i % 100000 == 0:
                    print(i, r1[:4])
            prev = r1[1]
    f.close()
        
    
if __name__ == "__main__":
    import proc_util
    import file_util
    fasta = "/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh38d1/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    # fasta = "/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.reorder.fasta"
    # fasta = "/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/hg38/Homo_sapiens_assembly38.fasta"
    # fasta = "/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/hg38.ucsc/ucsc.hg38.sorted.fa"
    # tsvgz = "/home/mk446/mutanno/DATASOURCE/PATHOGENICITY/CADD/hg38/v1.5/whole_genome_SNVs.tsv.gz"
    # baseidx = 2
    
    # tsvgz = "/home/mk446/mutanno/DATASOURCE/SPLICING/SpliceAI/hg38/spliceai_scores.raw.snv.hg38.vcf.gz"
    # baseidx = 3

    tsvgz = sys.argv[1]
    baseidx = int(sys.argv[2])
    check_base(tsvgz, baseidx)
