#!/usr/bin/env python
# -*- coding: utf-8 -*-
# make_datasource_based_on_vcf.py
# made by Daniel Minseok Kwon
# 2020-05-05 16:23:59
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


def make_datasource_based_on_vcf(vcf):

    k = 0
    for line in file_util.gzopen(vcf):
        line = file_util.decodeb(line)
        if line[0] != '#':
            k += 1
            arr = line.split('\t')
            chrom = arr[0]
            pos = int(arr[1])

            spos = pos
            epos = pos

            cmd = "python /home/mk446/bio/mutanno/SRC/src/mutanno.py makedata "
            cmd += "-ds /home/mk446/bio/mutanno/SRC/tests/data/datastructure_v0.4.4ds.json "
            cmd += "-out /home/mk446/bio/mutanno/DATASOURCE/MAINANNOT/tmp/mc_"+str(k)+".tsi "
            cmd += "-vartype SNV "
            cmd += "-blocksize 10 "
            cmd += "-region "+chrom+":"+str(spos)+"-"+str(epos)+" ;"
            # print(cmd)
            
            out = "/home/mk446/bio/mutanno/DATASOURCE/MAINANNOT/tmp/mc_"+str(k)+".tsi"

            out2 = "/home/mk446/bio/mutanno/DATASOURCE/MAINANNOT/mc_3k.tsi"
            if k == 1:
                cmd = "cat " + out + " > " + out2 + ";"
            else:
                cmd = "cat " + out + " | grep -v '^#' >> " + out2 + ";"
            print(cmd)
    

if __name__ == "__main__":
    import proc_util
    import file_util
    # vcf = "/home/mk446/mutanno/TEST/NOVO2_all_variants_jc50_wgenome.sorted.vcf.gz"
    vcf = "/home/mk446/mutanno/TEST/GAPFI7HPIVFM.mod.vcf"
    # vcf = "/home/mk446/mutanno/TEST/NOVO2_all_variants_jc50_wgenome.sorted.test.vcf"
    make_datasource_based_on_vcf(vcf)
