#!/usr/bin/env python
# -*- coding: utf-8 -*-
# split_run_datasource_4_microannot.py
# made by Min-Seok Kwon
# 2020-01-14 12:14:40
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
sys.path.append('..')

def print_log(log, cont):
    file_util.fileSave(log, cont + '\n', 'a')

# This function needs lots of computing time. 
def merge_and_tabixgz(path, out):
    log = out + '.log'
    file_util.fileSave(log, '', 'w')
    chromlist = seq_util.CHROM_LIST['b38d']
    
    # regionlist = seq_util.get_split_region()
    # print(regionlist)
    # print(len(regionlist))
    # for r1 in regionlist:

    for chrom in chromlist:
        if len(chrom) < 3:
            if chrom in ['1','2','3','4','5','6','7','8']:
                infile = path + 'all.tsi'
                if chrom == '1':
                    i = 0
                    for line in open(infile):
                        i += 1
                        if line[0] == '9':
                            break
                        else:
                            print(line, end = '')
                        if i % 10000000 == 0:
                            print_log(log, line[:20])
                            i = 0
                            # break
            else:
                infile = path + 'allchrom_' + chrom + '.tsi'
                i = 0
                for line in open(infile):
                    if (line[0] == '#' and chrom == '1') or line[0] != '#':
                        print(line, end = '')
                    i += 1
                    if i % 10000000 == 0:
                        print_log(log, line[:20])
                        i = 0
                        # break
            

def run(out):

    cmd = "python s04_merge_and_tabixgz.py merge | bgzip -c > " + out + ".gz; tabix -f -p vcf "+out+".gz"
    print(cmd)


if __name__ == "__main__":
    import proc_util
    import file_util
    import seq_util
    import preproc_util
    path = preproc_util.DATASOURCEPATH + "/MICROANNOT/tmp/"

    out = path + 'all2.tsi'
    if len(sys.argv) == 1:
        run(out)
    elif sys.argv[1] == 'merge':
        merge_and_tabixgz(path, out)

    # ds_json_file = "/home/mk446/mutanno/SRC/tests/data/datastructure_microannot_v0.4.1ds.json"
    # out = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/microannot_datasource.#CHROM#.v0.4.1_200421.tsi"
    # split_run_datasource_4_microannot(ds_json_file, path)
    
    # ds_json_file = "/home/mk446/mutanno/SRC/tests/data/datastructure_microannot_v0.4.1ds.json"
    # out = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/microannot_datasource.#CHROM#.v0.4.1_200421.tsi"
    # split_run_datasource_4_microannot_chrom(ds_json_file, out)
    
