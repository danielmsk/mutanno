#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s08_make_one_source.py
# made by Min-Seok Kwon
# 2020-01-21 06:23:00
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


def s08_make_one_source():
    for chrom in seq_util.MAIN_CHROM_LIST:
        if chrom == "MT":
            chrom = "M"
        for fname in file_util.walk(path + "chr" + chrom + "/", '.tsv.vcf.gz'):
            # print(fname)
            arr = fname.split('/')[-1].split('.')[0].split('_')
            newpath = '/'.join(fname.split('/')[:-1])
            cmd = "mutanno makedata -ds /home/mk446/mutanno/SRC/tests/datastructure_microannot_v1.0.json "
            cmd += "-out " + newpath + "/chr" + chrom + "_" + arr[1] + "_" + arr[2] + "_microannot "
            cmd += "-vartype SNP "
            cmd += "-region " + chrom + ":" + arr[1] + "-" + arr[2] + ""
            print(cmd)
            # break
        # break


def s08_merge_vep_by_chrom_pipe_bgzip(chrom):
    logfile = path + "microannot.hg38." + chrom + ".tsv.gz.log"
    file_util.fileSave(logfile, '', 'w')
    i = 0
    for k in range(400):
        vcfmap = {}
        for tsv in file_util.walk(path + "chr" + chrom + "/" + str(k) + "/", '_microannot.tsv'):
            k1 = int(tsv.split('/')[-1].split('_')[1])
            vcfmap[tsv] = k1
        (ks, vs) = struct_util.sortdict(vcfmap)

        for tsv in ks:
            j = 0
            log = str(i + 1) + ': ' + tsv
            file_util.fileSave(logfile, log + '\n', 'a')
            for line in file_util.gzopen(tsv):
                # line = line.decode('UTF-8')
                if (i == 0 and line[0] == '#'):
                    line = line.replace('\t.\t', '\tID\t')
                if (i == 0) or (i > 0 and line[0] != '#'):
                    print(line, end='')
                j += 1
                if j % 100000 == 0:
                    file_util.fileSave(logfile, '\t' + line[:30] + '\n', 'a')
            i += 1

    file_util.fileSave(logfile, 'Done.\n', 'a')


def run():
    for chrom in seq_util.MAIN_CHROM_LIST:
        if chrom == "MT":
            chrom = "M"
        outgz = path + "microannot_datasource." + chrom + ".v1.0.tsv.gz"
        cmd = "python /home/mk446/mutanno/SRC/scripts/precal_vep/s08_make_one_source.py " + chrom
        cmd += " | bgzip -c > " + outgz + ";"
        cmd += "sleep 120;"
        cmd += "tabix -f -p vcf " + outgz + ";"
        runsh = 's08_chrom' + chrom + '.sh'
        file_util.fileSave(runsh, cmd + '\n', 'w')
        print('Saved..', runsh)


if __name__ == "__main__":
    import proc_util
    import file_util
    import seq_util
    import struct_util
    path = "/home/mk446/mutanno/PRECALVEP/"
    spath = "/home/mk446/mutanno/SRC/scripts/precal_vep/"

    # s08_make_one_source()

    if len(sys.argv) == 1:
        run()
    else:
        s08_merge_vep_by_chrom_pipe_bgzip(sys.argv[1])
