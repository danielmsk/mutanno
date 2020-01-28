#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s01_convert_vcf2tsv.py
# made by Daniel Minseok Kwon
# 2020-01-13 16:13:05
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


def rm_empty_vcf():
    out = "s02_rm_empty_vcf.sh"
    f = open(out, 'w')
    for vcf in file_util.walk(path, '.vcf'):
        # print(done)
        done = vcf + '.vep.txt.done'
        vep = done.replace('.vcf.vep.txt.done', '.vcf.vep.txt')
        vephtml = done.replace('.vcf.vep.txt.done', '.vcf.vep.txt_summary.html')
        varcnt = 0
        if file_util.is_exist(vcf):
            for line in file_util.gzopen(vcf):
                # line = line.decode('UTF-8')
                if line[0] != "#":
                    varcnt += 1
                    if varcnt > 2:
                        break
            if varcnt == 0 and not file_util.is_exist(vep) and not file_util.is_exist(vephtml):
                cmd = "rm " + vcf + '*;\n'
                proc_util.run_cmd(cmd, True)
                # print(cmd)
                # file_util.fileSave(out, cmd, 'a')
                f.write(cmd)
        else:
            print("Error (no vcf file):", vcf)
    f.close()
    proc_util.run_cmd('chmod 755 ' + out)


def s02_1_rm_all_done():
    print('run s02_1_rm_all_done')
    for done in file_util.walk(path, '.vcf.vep.txt.done'):
        vcf = done.replace('.vcf.vep.txt.done', '.vcf')
        vep = done.replace('.vcf.vep.txt.done', '.vcf.vep.txt')
        tsv = done.replace('.vcf.vep.txt.done', '.tsv')
        tsvgz = done.replace('.vcf.vep.txt.done', '.tsv.gz')
        tsvgztbi = done.replace('.vcf.vep.txt.done', '.tsv.gz.tbi')
        tsvgzchecked = done.replace('.vcf.vep.txt.done', '.tsv.gz.checked')
        tsvdone = done.replace('.vcf.vep.txt.done', '.tsv.done')
        vephtml = done.replace('.vcf.vep.txt.done', '.vcf.vep.txt_summary.html')
        vepchecked = done.replace('.vcf.vep.txt.done', '.vcf.vep.txt.checked')

        if file_util.is_exist(vep) and file_util.is_exist(vephtml) and file_util.is_exist(vepchecked):
            if file_util.is_exist(tsv) and file_util.is_exist(tsvgz) and file_util.is_exist(tsvgztbi):
                if file_util.is_exist(tsvgzchecked):
                    cmd = "rm " + vcf + "*;"
                    cmd += "rm " + tsv + ";"
                    cmd += "rm " + tsvdone + ";"
                    cmd += "rm " + tsv + ".gz.checked;"
                    proc_util.run_cmd(cmd, True)


def s02_2_check_tsvgzchecked():
    out = path + "aas02_2_check_tsvgzchecked.sh"
    print('run s02_2_check_tsvgzchecked')
    if "aas02_2_" not in runjoblist:
        k = 0
        for done in file_util.walk(path, '.vcf.vep.txt.done'):
            vcf = done.replace('.vcf.vep.txt.done', '.vcf')
            vep = done.replace('.vcf.vep.txt.done', '.vcf.vep.txt')
            tsv = done.replace('.vcf.vep.txt.done', '.tsv')
            tsvgz = done.replace('.vcf.vep.txt.done', '.tsv.gz')
            tsvgztbi = done.replace('.vcf.vep.txt.done', '.tsv.gz.tbi')
            tsvgzchecked = done.replace('.vcf.vep.txt.done', '.tsv.gz.checked')
            vephtml = done.replace('.vcf.vep.txt.done', '.vcf.vep.txt_summary.html')
            vepchecked = done.replace('.vcf.vep.txt.done', '.vcf.vep.txt.checked')

            if file_util.is_exist(vep) and file_util.is_exist(vephtml) and file_util.is_exist(vepchecked):
                if file_util.is_exist(tsv) and file_util.is_exist(tsvgz) and file_util.is_exist(tsvgztbi):
                    if not file_util.is_exist(tsvgzchecked):
                        k += 1
                        cmd = "python /home/mk446/mutanno/SRC/mutanno.py precal -check_vep_result "
                        cmd += "-vcf " + vcf + " "
                        cmd += "-vep_result " + tsvgz + ";"
                        out2 = out + '_' + str(k) + '.sh'
                        file_util.fileSave(out2, cmd, 'w')
                        proc_util.run_cmd('chmod 755 ' + out2)
                        print('\t', out2)

        time.sleep(5)
        proc_util.run_cmd('mv ' + out + '_*.sh /home/mk446/jobs/.')


def s02_3_tabixgz_tsv():
    global runjoblist
    out = path + "aas02_3_tabixgz_tsv.sh"
    print('run s02_3_tabixgz_tsv')
    if "aas02_3_" not in runjoblist:
        k = 0
        for done in file_util.walk(path, '.vcf.vep.txt.done'):
            # vcf = done.replace('.vcf.vep.txt.done', '.vcf')
            vep = done.replace('.vcf.vep.txt.done', '.vcf.vep.txt')
            tsv = done.replace('.vcf.vep.txt.done', '.tsv')
            tsvdone = done.replace('.vcf.vep.txt.done', '.tsv.done')
            tsvgz = done.replace('.vcf.vep.txt.done', '.tsv.gz')
            tsvgztbi = done.replace('.vcf.vep.txt.done', '.tsv.gz.tbi')
            tsvgzchecked = done.replace('.vcf.vep.txt.done', '.tsv.gz.checked')
            vephtml = done.replace('.vcf.vep.txt.done', '.vcf.vep.txt_summary.html')
            vepchecked = done.replace('.vcf.vep.txt.done', '.vcf.vep.txt.checked')

            if file_util.is_exist(vep) and file_util.is_exist(vephtml) and file_util.is_exist(vepchecked):
                if file_util.is_exist(tsv) and not file_util.is_exist(tsvgz) and not file_util.is_exist(tsvgztbi):
                    if not file_util.is_exist(tsvgzchecked):
                        if file_util.is_exist(tsvdone):
                            k += 1
                            cmd = "tabixgz " + tsv
                            out2 = out + '_' + str(k) + '.sh'
                            file_util.fileSave(out2, cmd, 'w')
                            print('\t', out2)
                            proc_util.run_cmd('chmod 755 ' + out2)

        time.sleep(5)
        proc_util.run_cmd('mv ' + out + '_*.sh /home/mk446/jobs/.')


def s02_4_conv_vep2tsv():
    global runjoblist
    out = path + "aas02_4_conv_vep2tsv.sh"
    print('run s02_4_conv_vep2tsv')

    if "aas02_4_" not in runjoblist:
        k = 0
        for vepchecked in file_util.walk(path, '.vcf.vep.txt.checked'):
            # vcf = done.replace('.vcf.vep.txt.done', '.vcf')
            done = vepchecked.replace('.vcf.vep.txt.checked', '.vcf.vep.txt.done')
            vep = done.replace('.vcf.vep.txt.done', '.vcf.vep.txt')
            tsv = done.replace('.vcf.vep.txt.done', '.tsv')
            tsvgz = done.replace('.vcf.vep.txt.done', '.tsv.gz')
            tsvdone = done.replace('.vcf.vep.txt.done', '.tsv.done')
            tsvgztbi = done.replace('.vcf.vep.txt.done', '.tsv.gz.tbi')
            tsvgzchecked = done.replace('.vcf.vep.txt.done', '.tsv.gz.checked')
            vephtml = done.replace('.vcf.vep.txt.done', '.vcf.vep.txt_summary.html')
            # vepchecked = done.replace('.vcf.vep.txt.done', '.vcf.vep.txt.checked')

            if file_util.is_exist(vep) and file_util.is_exist(vepchecked):
                if not file_util.is_exist(tsv) and not file_util.is_exist(tsvgz) and not file_util.is_exist(tsvgztbi):
                    if not file_util.is_exist(tsvgzchecked):
                        k += 1
                        cmd = mutanno + "convert -vep2tab -in " + vep + " -out " + tsv + ";"
                        cmd += "touch " + tsvdone + ";"
                        out2 = out + '_' + str(k) + '.sh'
                        file_util.fileSave(out2, cmd, 'w')
                        print('\t', out2)
                        proc_util.run_cmd('chmod 755 ' + out2)

        time.sleep(5)
        proc_util.run_cmd('mv ' + out + '_*.sh /home/mk446/jobs/.')


def s02_5_check_vep():
    global runjoblist
    out = path + "aas02_5_check_vep.sh"
    # f = open(out, 'w')
    print('run s02_5_check_vep')
    if "aas02_5_" not in runjoblist:
        k = 0
        for done in file_util.walk(path, '.vcf.vep.txt.done'):
            vcf = done.replace('.vcf.vep.txt.done', '.vcf')
            vep = done.replace('.vcf.vep.txt.done', '.vcf.vep.txt')
            tsv = done.replace('.vcf.vep.txt.done', '.tsv')
            tsvgz = done.replace('.vcf.vep.txt.done', '.tsv.gz')
            tsvgztbi = done.replace('.vcf.vep.txt.done', '.tsv.gz.tbi')
            tsvgzchecked = done.replace('.vcf.vep.txt.done', '.tsv.gz.checked')
            vephtml = done.replace('.vcf.vep.txt.done', '.vcf.vep.txt_summary.html')
            vepchecked = done.replace('.vcf.vep.txt.done', '.vcf.vep.txt.checked')

            if file_util.is_exist(vep) and file_util.is_exist(vephtml) and not file_util.is_exist(vepchecked):
                if not file_util.is_exist(tsv) and not file_util.is_exist(tsvgz) and not file_util.is_exist(tsvgztbi):
                    if not file_util.is_exist(tsvgzchecked):
                        k += 1
                        cmd = "python /home/mk446/mutanno/SRC/mutanno.py precal -check_vep_result "
                        cmd += "-vcf " + vcf + " "
                        cmd += "-vep_result " + vep + ";"
                        # f.write(cmd + '\n')
                        out2 = out + '_' + str(k) + '.sh'
                        file_util.fileSave(out2, cmd, 'w')
                        print('\t', out2)
                        proc_util.run_cmd('chmod 755 ' + out2)
        # f.close()
        # proc_util.run_cmd('chmod 755 ' + out)
        time.sleep(5)
        proc_util.run_cmd('mv ' + out + '_*.sh /home/mk446/jobs/.')


def get_corelist():
    corelist = []
    cmd = 'squeue -a -u mk446 --format="%.18i %.9P %.8u %.8T %.10M %.9l %.6D %R %j" -t R'
    # cmd = 'squeue -a --format="%.18i %.9P %.8u %.8T %.10M %.9l %.6D %R %j" -t R'
    part_cont = proc_util.run_cmd(cmd)
    part_list = part_cont.strip().split('\n')
    for line in part_list:
        if 'gatk_hc_' in line:
            cid = line.split('gatk_hc_')[-1].strip()
            # print (cid)
            corelist.append(cid)
    return corelist


def get_runjoblist():
    corelist = get_corelist()
    runjoblist = ""
    for jid in corelist:
        k = 0
        for shfile in file_util.walk(JOBPATH + jid + '/', '.sh'):
            if not file_util.is_exist(shfile + '.done'):
                # print (jid, shfile , "RUNNING...")
                runjoblist += shfile + '\n'
                # print(log)
                k += 1
    for fname in file_util.listdir(JOBPATH, '.sh'):
        runjoblist += fname + '\n'

    return runjoblist


if __name__ == "__main__":
    import proc_util
    import file_util

    JOBPATH = '/home/mk446/jobs/'
    path = "/home/mk446/mutanno/PRECALVEP/"
    spath = "/home/mk446/mutanno/SRC/scripts/precal_vep/"
    mutanno = "python /home/mk446/mutanno/SRC/mutanno.py "
    # rm_empty_vcf()

    while True:
        runjoblist = get_runjoblist()

        # s02_1_rm_all_done()
        # s02_2_check_tsvgzchecked()
        # s02_3_tabixgz_tsv()
        # s02_4_conv_vep2tsv()
        s02_5_check_vep()
        # print('....')
        # time.sleep(120)
        break
