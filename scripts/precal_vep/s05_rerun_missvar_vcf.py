#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s05_rerun_missvar_vcf.py
# made by Min-Seok Kwon
# 2020-01-15 17:24:46
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


def s05_vep_split_run():
    for line in open(missvcflist):
        vcf = line.strip()
        vep = vcf + ".vep.txt"
        # cmd = "vepsplit run " + vcf + " 1000;"
        # cmd = "rm -rf " + vcf + "*"
        # cmd = "mv " + vcf + "_tmp/*sh /home/mk446/jobs/."
        # cmd = "vepsplit merge " + vcf + " 1000;"
        # cmd = "python /home/mk446/mutanno/SRC/mutanno.py precal -check_vep_result -vcf " + vcf + " -vep_result " + vep
        # cmd = "rm -rf " + vep + ".checked"
        tsv = vcf.replace('.vcf', '') + '.tsv'
        if file_util.is_exist(tsv + ".gz.checked"):
            # cmd = "python /home/mk446/mutanno/SRC/mutanno.py precal -check_vep_result -vcf " + vcf + " -vep_result " + tsv + '.gz'
            # cmd = "tabixgz " + tsv
            cmd = "rm " + tsv + ".gz.checked"
            print(cmd)
            # proc_util.run_cmd(cmd)
        '''
        if not file_util.is_exist(vep + ".checked"):
            flag = True
            # print(vep)
            # print(cmd)
            # proc_util.run_cmd(cmd)
            for vcf2 in file_util.walk(vcf + '_tmp', '.vcf'):
                vep2 = vcf2 + ".vep.txt"
                cmd = "python /home/mk446/mutanno/SRC/mutanno.py precal -check_vep_result -vcf " + vcf2
                cmd += " -vep_result " + vep2
                # print(cmd)
                if not file_util.is_exist(vep2 + ".checked"):
                    # 231 miss
                    print(cmd)
                    flag = False
                    # print(vcf2)
                    if not file_util.is_exist(vcf2 + ".sh"):
                        cmd = "/home/mk446/bin/vep -i " + vcf2
                        cmd += " -o " + vcf2 + ".vep.txt --hgvs "
                        cmd += "--fasta /n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh38d1/GRCh38_full_analysis_set_plus_decoy_hla.fa "
                        cmd += "--assembly GRCh38 --use_given_ref --offline --cache_version 98 "
                        cmd += "--dir_cache /home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/nonindexed_vep_cache/homo_sapiens_merged "
                        cmd += "--plugin MaxEntScan,/home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/VEP_plugins-release-98/MaxEntScan/fordownload "
                        cmd += "--plugin TSSDistance --everything --force_overwrite --tab;"
                        cmd += "touch " + vcf2 + ".vep.txt.done;\n"
                        # file_util.fileSave(vcf2 + ".sh", cmd, 'w')
                        # proc_util.run_cmd("chmod 755 " + vcf2 + ".sh", True)
                    # if "chr18_56000001_56100000.vcf" not in vcf2 and "chr18_56100001_56200000.vcf" not in vcf2:
                    # cmd = "cp " + vcf2 + ".sh /home/mk446/jobs/."
                    # proc_util.run_cmd(cmd, True)
            if flag:
                pass
                # print(vcf)
                # cmd = "vepsplit merge " + vcf + " 100;"
                # print(cmd)
        # print(proc_util.run_cmd(cmd))
        '''


if __name__ == "__main__":
    import proc_util
    import file_util
    missvcflist = "/home/mk446/mutanno/PRECALVEP/missvcf.txt"
    s05_vep_split_run()
