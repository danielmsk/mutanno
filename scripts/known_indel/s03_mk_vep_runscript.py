#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s03_mk_vep_runscript.py
# made by Min-Seok Kwon
# 2020-01-14 11:36:44
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


def s03_mk_vep_runscript():
    i = 0
    k = 0
    header = ''
    f = ''
    for line in file_util.gzopen(vcf):
        line = line.decode('UTF-8')
        if line[0] == '#':
            header = line
        else:
            if i % varsize == 0:
                try:
                    f.close()
                except AttributeError:
                    pass
                k += 1

                no = str_util.zero_format(k, 6)
                inputvcf = path + title + '_' + no + '.vcf'
                shcmd = path + title + '_' + no + '.vcf.sh'
                # vep = path + title + '_' + no + '.vcf.vep.txt'
                vep = path + title + '_' + no + '.vcf.vep.vcf'
                f = open(inputvcf, 'w')
                print(inputvcf)
                f.write(header)

                cmd = "/home/mk446/bin/vep -i " + inputvcf + " -o " + vep + " --hgvs "
                cmd += "--fasta " + fasta + " --assembly GRCh38 --use_given_ref "
                cmd += "--offline --cache_version 98 --dir_cache " + vepcache + " "
                cmd += "--plugin MaxEntScan,/home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/VEP_plugins-release-98/MaxEntScan/fordownload "
                cmd += "--plugin TSSDistance "
                # cmd += "--everything --force_overwrite --tab;\n"
                cmd += "--everything --force_overwrite --vcf;\n"

                file_util.fileSave(shcmd, cmd, 'w')
                proc_util.run_cmd('chmod 755 ' + shcmd)
            i += 1
            f.write(line)
    f.close()


if __name__ == "__main__":
    import proc_util
    import file_util
    import str_util
    title = 'aaakid'
    varsize = 20000
    fasta = "/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh38d1/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    # vepcache = "/home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/nonindexed_vep_cache/homo_sapiens_merged"
    vepcache = "/home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/nonindexed_vep_cache/homo_sapiens_vep"
    vcf = '/home/mk446/mutanno/DATASOURCE/KNOWN_INDEL/hg38/known_indel.sorted.uniq.vcf.gz'
    path = '/home/mk446/mutanno/DATASOURCE/KNOWN_INDEL/hg38/tmp3/'
    s03_mk_vep_runscript()
