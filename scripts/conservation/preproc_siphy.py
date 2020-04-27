#!/usr/bin/env python
# -*- coding: utf-8 -*-
# preproc_siphy.py
# made by Daniel Minseok Kwon
# 2020-04-20 23:39:46
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
from pyliftover import LiftOver
sys.path.append('..')

lo = LiftOver('hg19', 'hg38')

def preproc_siphy(hg19file, out):
    
    f = open(out, 'w')
    fl = open(out + '.log', 'w')
    i = 0
    for line in file_util.gzopen(hg19file):
        i += 1
        line = file_util.decodeb(line)
        arr = line.split('\t')
        chrom = arr[0].strip()
        spos = int(arr[1].strip())
        epos = int(arr[2].strip())
        lift_spos = lo.convert_coordinate(chrom, spos)
        lift_epos = lo.convert_coordinate(chrom, epos)
        len1 = epos - spos
        # print(len(lift_spos),len(lift_epos), arr, lift_spos, lift_epos)
        log = []
        for k in range(len(lift_spos)):
            log = []
            lschrom = lift_spos[k][0]
            lspos = lift_spos[k][1]
            
            try:
                lift_epos[k]
                lepos = lift_epos[k][1]
                lechrom = lift_epos[k][0]
            except IndexError:
                lepos = lspos + len1
                log.append("NA_lift_epos")
                lechrom = lschrom

            if len(lift_spos) > 1:
                log.append("multipos")
            # if (lepos - lspos) == (epos-spos) and lschrom == lechrom:
                
            if lschrom != lechrom:
                log.append("diff_chrom")
            if (lepos - lspos) != (epos-spos):
                log.append("diff_poslen")
            cont = [lschrom, str(lspos), str(lepos), arr[3].strip(), arr[4].strip(), arr[0], arr[1], arr[2], ','.join(log)]
            f.write('\t'.join(cont) + '\n')

            if len(log) >= 1:
                cont = [lschrom, str(lspos), str(lepos), arr[0], arr[1], arr[2], ','.join(log)]
                fl.write('\t'.join(cont) + '\n')

        if len(lift_spos) == 0:
            log.append('NA_lift')
            cont = ['', '', '', arr[0], arr[1], arr[2], ','.join(log)]
            fl.write('\t'.join(cont) + '\n')

        if i % 1000000 == 0:
            print(i, arr)
            # break
    f.close()
    fl.close()
    print('Saved',out)
    proc_util.run_cmd('vcf-sort -c ' + out + ' > ' + out + '.sorted.bed', True)

def liftover_test():
    print(lo.convert_coordinate('chr1', 12677863))

if __name__ == "__main__":
    import proc_util
    import file_util
    import preproc_util
    path = preproc_util.DATASOURCEPATH + "/CONSERVATION/SIPHY/29mammals/"
    pifile = path + "hg19_29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.gz"
    outpifile = path + "hg38liftover_29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.tsv"
    preproc_siphy(pifile, outpifile)
    # liftover_test()

    omegafile = path + "hg19_29way_omega_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.gz"
    outomegafile = path + "hg38liftover_29way_omega_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt"
    preproc_siphy(omegafile, outomegafile)
