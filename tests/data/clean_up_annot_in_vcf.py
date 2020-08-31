#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### clean_up_annot_in_vcf.py
#### made by Daniel Minseok Kwon
#### 2020-08-31 12:12:24
#########################
import sys
import os
import time
from mutanno.util import file_util
from mutanno.util import proc_util

BASIC_INFO_FIELDS = ['AC', 'AF', 'AN', 'BaseQRankSum', 'DP',
                     'ExcessHet', 'FS', 'MLEAC', 'MLEAF', 'MQ', 'MQRankSum', 'QD', 'ReadPosRankSum', 'SOR', 'RAW_MQandDP']

def clean_up_annot_in_vcf(vcf_file):
    out = vcf_file + '.cleanup.vcf'
    f = open(out, 'w')

    rejected = {}
    for line in file_util.gzopen(vcf_file):
        line = file_util.decodeb(line)
        if line[0] == "#":
            flag = True
            if "##MUTANNO=" in line:
                flag = False
            if "##INFO=<ID=" in line:
                f1 = line.replace('##INFO=<ID=', '').split(',')[0]
                if f1 not in BASIC_INFO_FIELDS:
                    flag = False
            
            if flag:
                f.write(line)
        else:
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()

            infofields = []
            for f1 in arr[7].split(';'):
                arrf = f1.split('=')
                if len(arrf) == 1:
                    infofields.append(f1)
                else:
                    if arrf[0] in BASIC_INFO_FIELDS:
                        infofields.append(f1)
                    else:
                        try:
                            rejected[arrf[0]]
                        except KeyError:
                            print(arrf[0])
                            rejected[arrf[0]] = 1

            arr[7] = ';'.join(infofields)

            line = '\t'.join(arr) + '\n'
            f.write(line)
    f.close()

    time.sleep(1)
    cmd = "tabixgz " + out
    proc_util.run_cmd(cmd, True)

    print('Saved', out + '.gz')




if __name__ == "__main__":
    print('#USAGE: python clean_up_annot_in_vcf.py [VCF]')
    vcf_file = sys.argv[1]
    clean_up_annot_in_vcf(vcf_file)
