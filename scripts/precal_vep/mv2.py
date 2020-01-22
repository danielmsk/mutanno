#!/usr/bin/env python
# -*- coding: utf-8 -*-
# mv.py
# made by Min-Seok Kwon
# 2019-12-29 17:37:50
#########################
import time
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


def run():
    for d1 in file_util.listdir('./'):
        if 'chr' in d1:
            for d2 in file_util.listdir('./' + d1):
                # print(d1 + '/' + d2)
                cmd = "python " + path + "mv.py " + d1 + '/' + d2
                print(cmd)
                # for fname in file_util.listdir('./' + d1 + '/' + d2, '.vep.sh'):
                #     f1 = path + d1 + '/' + d2 + '/' + fname
                #     # print(f1)
                #     cmd = "mv " + f1 + ' ' + f1.replace('.vep.sh', '.vep.txt')
                #     # print(cmd)
                #     proc_util.run_cmd(cmd)
                # break
            # break


def mv(d1):
    for fname in file_util.listdir(path + d1, '.vcf'):
        f1 = path + d1 + '/' + fname
        # print(f1)
        vcf = path + d1 + '/' + fname
        vep = path + d1 + '/' + fname + '.vep.txt'
        cmd = "mv " + f1 + ' ' + f1.replace('.vep.sh', '.vep.txt')
        cmd = "python /home/mk446/mutanno/SRC/mutanno.py precal -check_vep_result -vcf " + vcf + " -vep_result " + vep
        # print(cmd)
        print(proc_util.run_cmd(cmd, True))


def s11_check_undone(d1):
    vcflist = []
    for fname in file_util.listdir(path + d1, '.vcf'):
        vcf = path + d1 + '/' + fname
        if not file_util.is_exist(vcf + '.vep.txt.checked'):
            cnt = 0
            for line in open(vcf):
                if line[0] != '#':
                    cnt += 1
            if cnt > 0:
                vcflist.append(vcf)
    out = path + 'checked_' + d1.replace('/', '_')
    file_util.fileSave(out, '\n'.join(vcflist), 'w')


def s12_merge_checked():
    for fname in file_util.listdir('./'):
        if 'checked_chr' in fname:
            # print(fname)
            if file_util.getFileSize(fname) > 0:
                flist = file_util.fileOpen(fname).split('\n')
                cmd = ""
                for inputvcf in flist:
                    out = inputvcf + '.vep.txt'
                    cmd += "/home/mk446/bin/vep -i " + inputvcf + " -o " + out + " --hgvs "
                    cmd += "--fasta " + fasta + " --assembly GRCh38 --use_given_ref "
                    cmd += "--offline --cache_version 98 --dir_cache " + vepcache + " "
                    cmd += "--plugin MaxEntScan,/home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/VEP_plugins-release-98/MaxEntScan/fordownload "
                    cmd += "--plugin TSSDistance "
                    cmd += "--everything --force_overwrite --tab;\n"
                print(cmd.strip())


def s13_getvcflist(d1):
    vcflist = []
    for fname in file_util.listdir(path + d1, '.vcf'):
        vcf = d1 + '/' + fname
        vcflist.append(vcf)
    out = path + 'vcf_' + d1.replace('/', '_')
    file_util.fileSave(out, '\n'.join(vcflist) + '\n', 'w')


def s14_merge_vcflist():
    for chrom in seq_util.MAIN_CHROM_LIST:
        # print(chrom)
        if chrom == "MT":
            chrom = "M"
        for k in range(1000):
            listfile = "./vcflist/vcf_chr" + chrom + "_" + str(k)
            if file_util.is_exist(listfile):
                if k == 0:
                    cmd = "cat " + listfile + " > vcflist_chr" + chrom
                else:
                    cmd = "cat " + listfile + " >> vcflist_chr" + chrom
                print(cmd)

                proc_util.run_cmd(cmd)


def s15_check_rerun():
    for line in open('r.sh'):
        arr = line.split(' ')
        vcf = arr[2]
        vep = arr[4]
        if not file_util.is_exist(vep):
            # print(vcf)
            cntvar = 0
            header = ''
            for line in open(vcf):
                if line[0] != '#':
                    arr = line.split('\t')
                    if arr[3].strip() != '':
                        cntvar += 1
                else:
                    header = line
            # print(header)
            if cntvar == 0:
                # print(vcf)
                # file_util.fileSave(vep + '.checked', '', 'w')
                # file_util.fileSave(vcf, header, 'w')
                pass
        else:
            if file_util.is_exist(vep + '.error'):
                # print(vcf)
                cmd = "rm " + vep + '.error'
                # print(proc_util.run_cmd(cmd))

                cmd = "python /home/mk446/mutanno/SRC/mutanno.py precal -check_vep_result -vcf " + vcf + " -vep_result " + vep
                # print(cmd)
                # print(proc_util.run_cmd(cmd))

                cmd = "/home/mk446/bin/vep -i " + vcf + " -o " + vep + " --hgvs "
                cmd += "--fasta " + fasta + " --assembly GRCh38 --use_given_ref "
                cmd += "--offline --cache_version 98 --dir_cache " + vepcache + " "
                cmd += "--plugin MaxEntScan,/home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/VEP_plugins-release-98/MaxEntScan/fordownload "
                cmd += "--plugin TSSDistance "
                cmd += "--everything --force_overwrite --tab;"
                print(cmd)
                pass
            else:
                # print(vcf)
                pass


def s16_vcfgz():
    for line in file_util.gzopen("vcflist.gz"):
        line = line.decode('UTF-8')
        # cmd = "tabixgz " + path + line.strip()
        vcf = path + line.strip()
        vep = vcf + '.vep.txt'
        # if file_util.is_exist(vcf) and file_util.is_exist(vcf+'.gz') and file_util.is_exist(vcf+'.gz.tbi'):
        #     cmd = "rm " + vcf
        #     # cmd = "rm " + path + line.strip() + ".vep.sh_summary.html"
        #     print(cmd)
        #     proc_util.run_cmd(cmd)

        if file_util.is_exist(vep) and file_util.is_exist(vep + '.checked') and file_util.is_exist(vep + '.done'):
            cmd = "gz " + vep
            print(cmd)


def s17_vep2tab():
    for line in file_util.gzopen("vcflist.gz"):
        line = line.decode('UTF-8')
        vcf = path + line.strip()
        vep = vcf + '.vep.txt'
        tab = vcf + '.vep.tab'
        # if file_util.is_exist(vep):
        if True:
            cmd = "python /home/mk446/mutanno/SRC/mutanno.py convert -vep2tab"
            cmd += " -in " + vep + '.gz'
            cmd += " -out " + tab
            cmd += ";\n"
            # cmd += "tabixgz " + tab + ";"
            print(cmd)
        # break


def s18_gz(d1):
    for fname in file_util.listdir(path + d1, '.vcf.gz'):
        vcf = path + d1 + '/' + fname
        vep = vcf + '.vep.txt'
        # if file_util.is_exist(vcf):
        #     cmd = "tabixgz " + vcf
        #     proc_util.run_cmd(cmd)
        #     print(cmd)

        # vep = path + d1 + '/' + fname
        # if file_util.is_exist(vep) and file_util.is_exist(vep + '.checked') and file_util.is_exist(vep + '.done'):
        #     cmd = "gz " + vep
        #     proc_util.run_cmd(cmd)
        #     print(cmd)
        # tab = vep.replace('.vep.txt', '.vep.tab')
        tab = vcf + '.vep.tab'
        if file_util.is_exist(tab):
            cmd = "tabixgz " + tab
            proc_util.run_cmd(cmd)
            print(cmd)
    print('sleep 60')
    time.sleep(60)
    for fname in file_util.listdir(path + d1, '.vcf.gz'):
        vcf = path + d1 + '/' + fname
        # vep = path + d1 + '/' + fname
        vep = vcf + '.vep.txt'
        tab = vcf + '.vep.tab'
        if file_util.is_exist(vcf) and file_util.is_exist(vcf + '.gz') and file_util.is_exist(vcf + '.gz.tbi'):
            cmd = "rm " + vcf
            proc_util.run_cmd(cmd)
            print(cmd)
        # tab = vep.replace('.vep.txt', '.vep.tab')
        # if file_util.is_exist(vep) and file_util.is_exist(vep+'.gz') and file_util.is_exist(vep + '.checked') and file_util.is_exist(vep + '.done'):
        #     cmd = "rm " + vep
        #     proc_util.run_cmd(cmd)
        #     print(cmd)
        if file_util.is_exist(tab) and file_util.is_exist(tab + '.gz') and file_util.is_exist(tab + '.gz.tbi'):
            cmd = "rm " + tab
            proc_util.run_cmd(cmd)
            print(cmd)


def s19_vep2tab(d1):
    flag_sleep = False
    for fname in file_util.listdir(path + d1, '.vcf.gz'):
        vcf = path + d1 + "/" + fname.replace('.vcf.gz', '.vcf')
        vep = vcf + '.vep.txt'
        tab = vcf + '.vep.tab'
        # print(vep + ".gz")
        if file_util.is_exist(vep + ".gz") and not (file_util.is_exist(tab + ".gz") and file_util.is_exist(tab + ".gz.tbi")):
            cmd = "python /home/mk446/mutanno/SRC/mutanno.py convert -vep2tab"
            cmd += " -in " + vep + '.gz'
            cmd += " -out " + tab
            cmd += ";\n"
            # cmd += "tabixgz " + tab + ";"
            print(cmd)
            proc_util.run_cmd(cmd)
            flag_sleep = True
    if flag_sleep:
        print('sleep 60')
        time.sleep(60)
        flag_sleep = False

    for fname in file_util.listdir(path + d1, '.vcf.gz'):
        vcf = path + d1 + "/" + fname.replace('.vcf.gz', '.vcf')
        vep = vcf + '.vep.txt'
        tab = vcf + '.vep.tab'
        if file_util.is_exist(tab) and not (file_util.is_exist(tab + ".gz") and file_util.is_exist(tab + ".gz.tbi")):
            cmd = "tabixgz " + tab + ";"
            print(cmd)
            proc_util.run_cmd(cmd)
            flag_sleep = True
    if flag_sleep:
        print('sleep 60')
        time.sleep(60)

    for fname in file_util.listdir(path + d1, '.vcf.gz'):
        vcf = path + d1 + "/" + fname.replace('.vcf.gz', '.vcf')
        vep = vcf + '.vep.txt'
        tab = vcf + '.vep.tab'
        if file_util.is_exist(tab) and file_util.is_exist(tab + ".gz") and file_util.is_exist(tab + ".gz.tbi"):
            cmd = "rm " + tab + ''
            print(cmd)
            proc_util.run_cmd(cmd)


def s20_check_veptabgz(d1):
    for fname in file_util.listdir(path + d1, '.vcf.gz'):
        vcf = path + d1 + "/" + fname.replace('.vcf.gz', '.vcf')
        tab = vcf + '.vep.tab'
        if file_util.is_exist(tab + ".gz") and file_util.is_exist(tab + ".gz.tbi"):
            cmd = "python /home/mk446/mutanno/SRC/mutanno.py precal -check_vep_result "
            cmd += "-vcf " + vcf + ".gz "
            cmd += "-vep_result " + tab + ".gz "
            print(cmd)
            print(proc_util.run_cmd(cmd))
            # break
        else:
            cntvar = 0
            for line in file_util.gzopen(vcf + '.gz'):
                line = line.decode('UTF-8')
                if line[0] != '#':
                    cntvar += 1
            if cntvar == 0:
                cmd = "rm " + vcf + "*"
                print(cmd)
                proc_util.run_cmd(cmd)
            else:
                pass
                # print(cmd)
                # print(tab)


def s21_check_veptabgz2(d1):
    out = path + d1.replace('/', '_') + ".log.sh"
    # file_util.fileSave(out, '', 'w')
    for fname in file_util.listdir(path + d1, '.vcf.gz'):
        vcf = path + d1 + "/" + fname.replace('.vcf.gz', '.vcf')
        tab = vcf + '.vep.tab'
        vep = vcf + '.vep.txt'
        if file_util.is_exist(tab + ".gz") and file_util.is_exist(tab + ".gz.tbi") and file_util.is_exist(tab + ".gz.checked"):
            pass
        else:
            cmd = "rm " + vcf + ".vep*;"

            cmd += "/home/mk446/bin/vep -i " + vcf + ".gz -o " + vep + " --hgvs "
            cmd += "--fasta " + fasta + " --assembly GRCh38 --use_given_ref "
            cmd += "--offline --cache_version 98 --dir_cache " + vepcache + " "
            cmd += "--plugin MaxEntScan,/home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/VEP_plugins-release-98/MaxEntScan/fordownload "
            cmd += "--plugin TSSDistance "
            cmd += "--everything --force_overwrite --tab;"
            cmd += "sleep 5;"

            cmd += mutanno + "convert -vep2tab -in " + vep + " -out " + tab + ";"
            # proc_util.run_cmd(cmd, True)
            # file_util.fileSave(out, vcf + '\n', 'a')
            cmd += "sleep 5;"
            cmd += "tabixgz " + tab + ";"
            cmd += "sleep 5;"
            cmd += mutanno + "precal -check_vep_result -vep_result " + tab + ".gz -vcf " + vcf + ".gz;"
            cmd += "sleep 5;"
            file_util.fileSave(out, cmd + '\n', 'a')


def s22_rm_emptyvcf(d1):
    out = path + d1.replace('/', '_') + '.log.sh'
    for fname in file_util.listdir(path + d1, '.vcf.gz'):
        cnt = 0
        for line in file_util.gzopen(path + d1 + '/' + fname):
            line = line.decode('UTF-8')
            if line[0] != "#":
                cnt += 1
                if cnt > 2:
                    break
        if cnt == 0:
            cmd = "rm " + path + d1 + '/' + fname[:-3] + '*;\n'
            # print(cmd)
            file_util.fileSave(out, cmd, 'a')


def s23_merge_chrom(chrom):
    out = path + "vep.hg38." + chrom + ".tsv"
    # f = open(out, 'w')
    i = 0
    for k in range(300):
        vcfmap = {}
        for vcf in file_util.walk(path + "chr" + chrom + "/" + str(k) + "/", '.vcf.gz'):
            # print(vcf)
            k1 = int(vcf.split('/')[-1].split('_')[1])
            vcfmap[vcf] = k1
        # print(vcfmap)
        (ks, vs) = struct_util.sortdict(vcfmap)
        # print(ks)
        # print(vs)
        for vcf in ks:
            vep = vcf[:-3] + ".vep.tab.gz"
            if i == 0:
                cmd = "zcat " + vep + " > " + out
            else:
                cmd = "zcat " + vep + " | grep -v '^#' >> " + out
            print(cmd)
            i += 1
        # break
    # f.close()
    cmd = "sleep 20;"
    print(cmd)
    cmd = "tabixgz " + out
    print(cmd)


def run_chrom():
    for chrom in seq_util.MAIN_CHROM_LIST:
        if chrom == "MT":
            chrom = "M"
            out = path + "merge_" + chrom + ".sh"
            cmd = "python " + path + "mv.py " + chrom + " > " + out + ";"
            cmd += "sleep 20;"
            cmd += "mv "+out+" /home/mk446/jobs/.;"
            print(cmd)


def s24_rm_intermediate_files(d1):
    out = path + d1.replace('/', '_') + '.log.sh'
    for fname in file_util.listdir(path + d1, '.vcf.gz'):
        vcf = path + d1 + '/' + fname[:-3]
        cmd = "rm " + vcf + ".vep.sh_summary.html;"
        cmd += "rm " + vcf + ".vep.tab.gz.checked;"
        cmd += "rm " + vcf + ".vep.txt.checked;"
        cmd += "rm " + vcf + ".vep.txt.done;"
        cmd += "rm " + vcf + ".vep.txt.gz;"
        proc_util.run_cmd(cmd, True)


if __name__ == "__main__":
    import proc_util
    import file_util
    import seq_util
    import struct_util
    fasta = "/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh38d1/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    vepcache = "/home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/nonindexed_vep_cache/homo_sapiens_merged"
    path = "/home/mk446/mutanno/PRECALVEP/"
    mutanno = "python /home/mk446/mutanno/SRC/mutanno.py "
    # run()
    # mv(sys.argv[1])
    # s11_check_undone(sys.argv[1])
    # s12_merge_checked()
    # s13_getvcflist(sys.argv[1])
    # s14_merge_vcflist()
    # s15_check_rerun()
    # s16_vcfgz()
    # s17_vep2tab()
    # s18_gz(sys.argv[1])
    # s19_vep2tab(sys.argv[1])
    # s20_check_veptabgz(sys.argv[1])
    # s21_check_veptabgz2(sys.argv[1])
    # s22_rm_emptyvcf(sys.argv[1])
    s23_merge_chrom(sys.argv[1])
    # run_chrom()
    # s24_rm_intermediate_files(sys.argv[1])
