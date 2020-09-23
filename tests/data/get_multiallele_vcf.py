#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### get_multiallele_vcf.py
#### made by Daniel Minseok Kwon
#### 2020-08-09 13:34:51
#########################
import sys
import os
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path="/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path="/ms1/bin/python_lib"
else:
    sys_path="/home/mk446/bin/python_lib"
sys.path.append(sys_path)

import file_util
import proc_util


def get_multiallele_vcf(vcf, out):
    f = open(out, 'w')
    for line in file_util.gzopen(vcf):
        line = file_util.decodeb(line)
        if line[0] == "#":
            f.write(line)
        else:
            arr = line.split('\t')
            if ',' in arr[4]:
                f.write(line)
    f.close()
    proc_util.run_cmd('tabixgz ' + out, True)


if __name__ == "__main__":
    # vcf = "test_trio_1000_1.vcf.gz"
    # out = "test_trio_multiallele.vcf"
    vcf = "vep_output2.vcf.bak"
    out = "vep_output2_multiallele.vcf"
    get_multiallele_vcf(vcf, out)
