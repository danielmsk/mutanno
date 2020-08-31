#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### mk_target_sourcefile.py
#### made by Daniel Minseok Kwon
#### 2020-08-31 10:34:01
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
import tabix


def mk_target_sourcefile(ori_source_file, vcf_file):
    out = ori_source_file + '.target.mti'
    f = open(out, 'w')

    # save header
    for line in file_util.gzopen(ori_source_file):
        line = file_util.decodeb(line)
        if line[0] == '#':
            f.write(line)
            break

    m = {}
    tb = tabix.open(ori_source_file)
    for line in file_util.gzopen(vcf_file):
        line = file_util.decodeb(line)
        if line[0] != '#':
            arr = line.split('\t')
            chrom = arr[0].replace('chr', '')
            pos = int(arr[1])
            ref = arr[3]
            alt = arr[4]
            spos = pos
            epos = pos + len(alt)
            varkey = '_'.join(arr[:5])
            print(varkey)

            recs = tb.querys(chrom + ":"+str(spos)+"-" + str(epos))
            for r1 in recs:
                try:
                    m[varkey]
                except KeyError:
                    f.write('\t'.join(r1))
                    m[varkey] = 1

    f.close()

    cmd = "vcf-sort -c " + out + ' > ' + out + '.sorted.mti'
    proc_util.run_cmd(cmd, True)
    cmd = "tabixgz " + out + '.sorted.mti'
    proc_util.run_cmd(cmd, True)

    print('Final source file:', out + '.sorted.mti.gz')



    


if __name__ == "__main__":
    print('#USAGE: python mk_target_sorucefile.py [Origianl source file] [VCF]')
    ori_source_file = sys.argv[1]
    vcf_file = sys.argv[2]
    mk_target_sourcefile(ori_source_file, vcf_file)
