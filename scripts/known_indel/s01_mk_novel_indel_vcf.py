#!/usr/bin/env python
# -*- coding: utf-8 -*-
# mk_known_indel_vcf.py
# made by Min-Seok Kwon
# 2020-01-13 19:32:02
#########################

import sys
import os
import tabix
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path = "/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path = "/ms1/bin/python_lib"
else:
    sys_path = "/home/mk446/bin/python_lib"
sys.path.append(sys_path)


def is_var_in_vep(veptb, chrom, pos, ref, alt):
    flag = False
    recs = veptb[chrom].querys(chrom + ":" + pos + "-" + pos)
    for r1 in recs:
        if r1[1] == pos and r1[3] == ref and r1[4] == alt:
            flag = True
            break
    return flag

def mk_known_indel_vcf(chrom, spos, epos, f, vcftb, vep, chromidx, posidx, refidx, altidx):
    # print("processing...", vcf)

    recs1 = vcftb.querys(chrom + ":" + str(spos) + "-" + str(epos))
    for arr in recs1:
        chrom = arr[chromidx].replace('chr', '')
        pos = arr[posidx]
        ref = arr[refidx].strip()
        alt = arr[altidx].strip()

        if not (len(ref) == 1 and len(alt) == 1):
            try:
                tmp = veptb[chrom]
            except KeyError:
                veptb[chrom] = tabix.open(vep.replace('#CHROM#',chrom))

            flag = is_var_in_vep(veptb, chrom, pos, ref, alt)
            if not flag:
                cont = 'chr' + arr[chromidx].replace('chr', '') + '\t' + arr[posidx] + \
                    '\t.\t' + arr[refidx] + '\t' + arr[altidx]
                cont += addfield
                f.write(cont + '\n')
        if ',' in arr[refidx] or ',' in arr[altidx]:
            print('Error:', arr, vcf)



def run():
    chunklist = seq_util.get_split_region(2000000, 'b38d')
    for c1 in chunklist:
        chrom = c1[0]
        spos = c1[1]
        epos = c1[2]
        k = c1[3]

        cmd = "python /home/mk446/bio/mutanno/SRC/scripts/known_indel/s01_mk_novel_indel_vcf.py"
        cmd += " " + chrom
        cmd += " " + str(spos)
        cmd += " " + str(epos)
        cmd += " " + str(k)
        print(cmd)


def merge():
    chunklist = seq_util.get_split_region(2000000, 'b38d')
    for c1 in chunklist:
        k = str(c1[3])
        out2 = out + "_vcf/" + k + ".vcf"
        if k == 1:
            cmd = "cat " + out2 + " > " + out
        else:
            cmd = "cat " + out2 + " | grep -v '^#' >> " + out
        print(cmd)



if __name__ == "__main__":
    import file_util
    import seq_util
    addfield = '\t.\t.\t.\tGT:AD:DP\t0/1:0,30:60'
    header = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1'
    path = "/home/mk446/bio/mutanno/DATASOURCE/KNOWN_INDEL/hg38/"
    vep = path + "known_indel.chr#CHROM#.vep.tsi.gz"
    out = path + "known_indel_new.vcf"
    vcf = "/home/mk446/bio/mutanno/DATASOURCE/dbSNP/hg38/GRCh38_latest_dbSNP_all.convchrom.onlyid.indel.vcf.gz"
    veptb = {}
    vcftb = tabix.open(vcf)
    
    if len(sys.argv) == 1:
        # run()
        merge()
    else:
        k = sys.argv[4]
        out2 = out + "_vcf/" + k + ".vcf"
        f = open(out2 , 'w')
        f.write(header + '\n')
        mk_known_indel_vcf(sys.argv[1],int(sys.argv[2]),int(sys.argv[3]), f, vcftb, vep, 0, 1, 3, 4)
        f.close()
        print("Saved",out2)
