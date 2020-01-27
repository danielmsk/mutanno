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


VEP_OPTION_CMD = "#VEP# -i #INPUTVCF# -o #OUT#"
VEP_OPTION_CMD += " --hgvs"
VEP_OPTION_CMD += " --fasta #FASTA#"
VEP_OPTION_CMD += " --assembly GRCh38"
VEP_OPTION_CMD += " --use_given_ref"
VEP_OPTION_CMD += " --offline"
VEP_OPTION_CMD += " --cache_version #CACHE_VERSION#"
VEP_OPTION_CMD += " --dir_cache #VEPCACHE#"
VEP_OPTION_CMD += " --everything"
VEP_OPTION_CMD += " --force_overwrite"
VEP_OPTION_CMD += " --vcf"
# PLUGIN
# VEP_OPTION_CMD += "--plugin LoF --plugin LoF,human_ancestor_fa:/home/mk446/BiO/Data/vep/human_ancestor.fa.gz "
VEP_OPTION_CMD += " --plugin MaxEntScan,/home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/VEP_plugins-release-99/fordownload"
VEP_OPTION_CMD += " --plugin TSSDistance"
VEP_OPTION_CMD += " --dir_plugins /home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/VEP_plugins-release-99"
# v99 added
VEP_OPTION_CMD += " --plugin SpliceRegion,Extended"
# VEP_OPTION_CMD += " --plugin LoF,loftee_path:/home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/VEP_plugins-release-99/loftee"


def s05_vep_split_run():
    for line in open(missvcflist):
        vcf = line.strip()
        out = vcf + ".vep.txt"

        vepcmd = VEP_OPTION_CMD
        vepcmd = vepcmd.replace('#VEP#', vep)
        vepcmd = vepcmd.replace('#INPUTVCF#', vcf)
        vepcmd = vepcmd.replace('#OUT#', out)
        vepcmd = vepcmd.replace('#FASTA#', fasta)
        vepcmd = vepcmd.replace('#VEPCACHE#', vepcache)
        vepcmd = vepcmd.replace('#CACHE_VERSION#', cache_version)
        cmd = vepcmd + ";"
        cmd += "sleep 60;"
        cmd += "tabixgz " + out + ";"
        cmd += "touch " + out + ".done;"

        cmd = "mutanno precal -check_vep_result "
        cmd += "-vcf " + vcf + " "
        cmd += "-vep_result " + out + ".gz;"
        print(cmd)

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
    vep = "/home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/ensembl-vep-release-99/vep"
    fasta = "/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh38d1/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    vepcache = "/home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/nonindexed_vep_cache/homo_sapiens_vep"
    cache_version = "99"
    missvcflist = "/home/mk446/mutanno/PRECALVEP/missvcf.txt"
    s05_vep_split_run()
