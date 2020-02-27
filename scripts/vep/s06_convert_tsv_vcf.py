#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s06_convert_tsv_vcf.py
# made by Min-Seok Kwon
# 2020-01-16 17:17:08
#########################

import sys
import os
import time
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path = "/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path = "/ms1/bin/python_lib"
else:
    sys_path = "/home/mk446/bin/python_lib"
sys.path.append(sys_path)


def s06_convert_tsv_vcf(tsv):
    vcf = tsv.replace('.tsv.gz', '.tsv.vcf')
    f = open(vcf, 'w')
    i = 0
    for line in file_util.gzopen(tsv):
        line = line.decode('UTF-8')
        arr = line.split('\t')
        i += 1
        if line[0] == '#':
            f.write('\t'.join(arr[:2]) + '\tID\t' + '\t'.join(arr[2:]))
        else:
            f.write('\t'.join(arr[:2]) + '\t\t' + '\t'.join(arr[2:]))
    f.close()

    time.sleep(3)
    proc_util.run_cmd('tabixgz ' + vcf)
    time.sleep(3)

    j = 0
    for line in file_util.gzopen(vcf + '.gz'):
        j += 1

    if i == j:
        proc_util.run_cmd('rm -rf ' + vcf)


def run():
    for tsv in file_util.walk(path, '.tsv.gz'):
        vcf = tsv.replace('.tsv.gz', '.tsv.vcf.gz.tbi')
        if not file_util.is_exist(vcf):
            cmd = "python /home/mk446/mutanno/SRC/scripts/precal_vep/s06_convert_tsv_vcf.py " + tsv
            print(cmd)


if __name__ == "__main__":
    import proc_util
    import file_util
    path = "/home/mk446/mutanno/PRECALVEP/"
    if len(sys.argv) == 1:
        run()
    else:
        s06_convert_tsv_vcf(sys.argv[1])
