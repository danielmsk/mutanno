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


def s03_mk_vep_runscript():
    i = 0
    k = 0
    header = ''
    f = ''
    for line in file_util.gzopen(vcf):
        if vcf.endswith('.gz'):
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

                spath = path + str(int(k / 1000)) + '/'

                inputvcf = spath + title + '_' + no + '.vcf'
                shcmd = spath + title + '_' + no + '.vcf.sh'
                vep = spath + title + '_' + no + '.vcf.vep.txt'
                # vep = spath + title + '_' + no + '.vcf.vep.vcf'
                file_util.check_dir(inputvcf)

                f = open(inputvcf, 'w')
                print(inputvcf)
                f.write(header)

                if not file_util.is_exist(vep + "_summary.html"):
                    cmd = VEP_OPTION_CMD
                    cmd = cmd.replace('#VEP#', "/home/mk446/bin/vep")
                    cmd = cmd.replace('#INPUTVCF#', inputvcf)
                    cmd = cmd.replace('#OUT#', vep)
                    cmd = cmd.replace('#FASTA#', fasta)
                    cmd = cmd.replace('#VEPCACHE#', vepcache)
                    cmd = cmd.replace('#CACHE_VERSION#', '99')
                    cmd = cmd + ";"
                    cmd += "touch " + inputvcf + ".vep.txt.done;"

                    file_util.fileSave(shcmd, cmd, 'w')
                    proc_util.run_cmd('chmod 755 ' + shcmd)
            i += 1
            f.write(line)
    f.close()


if __name__ == "__main__":
    import proc_util
    import file_util
    import str_util
    title = 'sindel'
    varsize = 5000
    fasta = "/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh38d1/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    # vepcache = "/home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/nonindexed_vep_cache/homo_sapiens_merged"
    vepcache = "/home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/nonindexed_vep_cache/homo_sapiens_vep"
    # vcf = '/home/mk446/bio/mutanno/DATASOURCE/KNOWN_INDEL/hg38/known_indel.sorted.uniq.vcf.gz'
    # vcf = '/home/mk446/bio/mutanno/DATASOURCE/KNOWN_INDEL/hg38/known_indel_new.vcf'
    vcf = '/home/mk446/bio/mutanno/DATASOURCE/KNOWN_INDEL/hg38/known_indel_new.uniq.vcf'
    path = '/home/mk446/bio/mutanno/DATASOURCE/KNOWN_INDEL/hg38/tmp/'
    s03_mk_vep_runscript()
