#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### s03_listover.py
#### made by Min-Seok Kwon
#### 2019-11-05 12:01:55
#########################
import sys
import os
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path="/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path="/ms1/bin/python_lib"
else:
    sys_path="/home/mk446/bin/python_lib"
sys.path.append(sys_path)
import file_util
import proc_util
from pyliftover import LiftOver
import tabix
import seq_util

def make_listover_map(b37, b38, tchrom, tspos, tepos, k1):
    lo = LiftOver('hg19', 'hg38')
    # lo = LiftOver('hg38', 'hg19')

    tb_b38 = tabix.open(b38)
    tb_b37 = tabix.open(b37)

    out = "/home/mk446/bio/mutanno/s03_liftover_map/"+tchrom+"/"+tchrom+"_"+str(k1)+".tsv"
    # out = "/home/mk446/bio/mutanno/s03_liftover_map/"+tchrom+"_"+str(k1)+".tsv"
    file_util.check_dir(out)
    f = open(out, 'w')
    cont = ['#CHROM','POS','REF','LIFTOVER_hg38']
    f.write('\t'.join(cont)+'\n')
    prev = ""
    i = 0

    recs = tb_b37.query(tchrom, tspos, tepos)
    for arr in recs:
        # print (arr)
        if arr[0] + ':' + arr[1] != prev:
            liftpos = lo.convert_coordinate('chr'+arr[0], int(arr[1]))
            # print (p2[0][0], p2[0][1], p2)
            # print (liftpos, arr)
            ref_b37 = arr[2]
            cont = [arr[0]]
            cont.append(arr[1])
            cont.append(ref_b37)

            arr2 = []
            for p2 in liftpos:
                recs = tb_b38.query(p2[0].replace('chr',''), int(p2[1]), int(p2[1]))
                for r1 in recs:
                    ref_b38 = r1[2]
                    break
                arr2.append(p2[0].replace('chr','') + '_'+str(p2[1])+"_"+ref_b38)
            cont.append(';'.join(arr2))
            f.write('\t'.join(cont)+'\n')
            i += 1
            if ref_b37 != ref_b38:
                print ("#ERRROR:",cont)
            if i % 10000 == 0:
                print (i, cont)
                i = 0
            if len(liftpos) != 1:
                print (arr, liftpos)
        prev = arr[0] + ':' + arr[1]


def run():
    sublist = seq_util.get_split_region(20000, 'b38d')
    # print (len(sublist))
    for r1 in sublist:
        # if r1[0] != "1":
        # c1 = int(r1[0])
        # if c1 > 1 and c1 < 6:
        # if c1 >= 6 and c1 < 13:
        if r1[0] == "X" or r1[0] == "Y" or r1[0] == "M" or int(r1[0]) >= 13:
            cmd = "python /home/mk446/bio/mutanno/s03_liftover.py " + r1[0] + " " + str(r1[1]) + " " + str(r1[2]) + " " + str(r1[3])
            print (cmd)


def check_files():
    prev = ""
    sublist = seq_util.get_split_region(20000, 'b38d')
    # print (len(sublist))
    for r1 in sublist:

        cmd = "python /home/mk446/bio/mutanno/s03_liftover.py " + r1[0] + " " + str(r1[1]) + " " + str(r1[2]) + " " + str(r1[3])
        tchrom = r1[0]
        k1 = r1[3]
        out = "/home/mk446/bio/mutanno/s03_liftover_map/"+tchrom+"/"+tchrom+"_"+str(k1)+".tsv"

        if prev != tchrom:
            print ('#CHROM', tchrom)

        if not file_util.is_exist(out):
            print (cmd)

        prev = tchrom

        # if tchrom == "4":
        #     break



if __name__ == "__main__":
    # vcf = sys.argv[1]
    b37 = "/home/mk446/BiO/Data/CADD/v1.4/whole_genome_SNVs.tsv.gz"
    b38 = "/home/mk446/BiO/Data/CADD/v1.4/GRCh38/whole_genome_SNVs.tsv.gz"
    if len(sys.argv) > 2:
        tchrom = sys.argv[1]
        tspos = int(sys.argv[2])
        tepos = int(sys.argv[3])
        k1 = int(sys.argv[4])
        make_listover_map(b37, b38, tchrom, tspos, tepos, k1)
    else:
        # run()
        check_files()
