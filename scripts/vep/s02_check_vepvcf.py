#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s02_check_vepvcf.py
# made by Daniel Minseok Kwon
# 2020-01-27 10:10:20
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


def s02_check_vepvcf():
    cnt = {}
    cnt['vcf'] = 0
    cnt['vep'] = 0
    cnt['vepgz'] = 0
    cnt['vepgztbi'] = 0
    cnt['vepgzchecked'] = 0
    cnt['no_variant_vcf'] = 0
    cnt['need_to_vep'] = 0
    for vcf in file_util.walk(path, '.vcf'):
        vep = vcf + '.vep.txt'
        vepgz = vcf + '.vep.txt.gz'
        vepgzchecked = vcf + '.vep.txt.gz.checked'
        vepgztbi = vcf + '.vep.txt.gz.tbi'
        if not file_util.is_exist(vepgzchecked) and file_util.is_exist(vepgztbi):
            cmd = mutanno + " precal -check_vep_result"
            cmd += " -vcf " + vcf
            cmd += " -vep_result " + vepgz + ";"
            print(cmd)
        if not file_util.is_exist(vepgzchecked) and not file_util.is_exist(vepgztbi) and file_util.is_exist(vep):
            cmd = "tabixgz " + vep + ';'
            print(cmd)
        if not file_util.is_exist(vepgzchecked) and not file_util.is_exist(vepgztbi) and not file_util.is_exist(vep):
            cnt_line = 0
            for line in open(vcf):
                if line[0] != '#':
                    cnt_line += 1
            if cnt_line == 0:
                cnt['no_variant_vcf'] += 1
            else:
                cnt['need_to_vep'] += 1
        if file_util.is_exist(vepgzchecked):
            cnt['vepgzchecked'] += 1
        if file_util.is_exist(vep):
            cnt['vep'] += 1
        if file_util.is_exist(vepgz):
            cnt['vepgz'] += 1
        if file_util.is_exist(vepgztbi):
            cnt['vepgztbi'] += 1

        if file_util.is_exist(vepgzchecked) and file_util.is_exist(vepgztbi) and file_util.is_exist(vep):
            cmd = "rm " + vep
            print(cmd)
        cnt['vcf'] += 1
    for k1 in cnt.keys():
        print('#', k1, cnt[k1])


if __name__ == "__main__":
    import proc_util
    import file_util
    path = "/home/mk446/mutanno/PRECALVEP/"
    # path = "/home/mk446/mutanno/PRECALVEP/chr22/"
    mutanno = "python /home/mk446/mutanno/SRC/mutanno.py "
    s02_check_vepvcf()
