#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s02_mk_check_tsi.py
# made by Daniel Minseok Kwon
# 2020-02-10 09:04:55
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
import tabix


def get_covered_region(chrom, t_spos, t_epos):
    cov_list = []
    spos = -9
    ipos = -9
    tb = tabix.open(prev_tsi.replace('#CHROM#', chrom))
    for rec in tb.query(chrom, t_spos, t_epos):
        pos = int(rec[1])
        if ipos == pos:
            pass
        elif ipos + 1 == pos:
            ipos += 1
        elif ipos + 1 < pos and ipos > 0:
            epos = ipos
            cov_list.append([spos, epos])
            
            spos = pos -1 
            ipos = pos

        if spos < 0:
            spos = pos - 1
            ipos = pos

    if spos + 1 < ipos:
        epos = ipos
        cov_list.append([spos, epos])

    return cov_list

def s02_mk_check_tsi(tsifile, region):
    chrom = region.split(':')[0]
    spos = int(region.split(':')[1].split('-')[0])
    epos = int(region.split(':')[1].split('-')[1])

    cov_list = get_covered_region(chrom, spos, epos)
    print(cov_list[:5])

    flag = True
    tpos = spos
    pos_flag = False
    for line in open(tsifile):
        if line[0] != '#':
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            print(arr)
            break
            # pos = int(arr[1])
            # if pos == tpos + 1 and pos_flag:
            #     pos_flag = True
            #     tpos = pos + 1
            # elif pos == tpos:
            #     pos_flag = True
            # else:
            #     pos_flag = False
            #     print("ERROR ", arr[0], arr[1], tpos)
            #     break
            # if not pos_flag:
            #     flag = False

    if flag:
        donefile = tsifile + ".done"
        print(donefile)
        # file_util.fileSave(donefile, '', 'w')


if __name__ == "__main__":
    import proc_util
    import file_util
    import seq_util
    
    # covered_file = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/v0.4.1_200429_chrom/microannot_datasource.#CHROM#.v0.4.1_200429.tsi.gz.covered.bed.gz"
    prev_tsi = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/v0.4.1_200429_chrom/microannot_datasource.#CHROM#.v0.4.1_200429.tsi.gz"
    
    
    tsifile = sys.argv[1]
    region = sys.argv[2]
    
    s02_mk_check_tsi(tsifile, region)
