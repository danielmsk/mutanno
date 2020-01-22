#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### liftover.py
#### made by Min-Seok Kwon
#### 2019-11-05 17:30:11
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

def liftover(tsv):

    tb_b38 = tabix.open(b38)


    lo = LiftOver('hg19', 'hg38')

    out = tsv.replace('.tsv.gz','') + '.liftover.tsv'
    out2 = tsv.replace('.tsv.gz','') + '.liftover.unmatch.tsv'
    i = 0
    f = open(out, 'w')
    f2 = open(out2, 'w')
    for line in file_util.gzopen(tsv):
        line = line.decode('UTF-8')
        if line[0] == "#":
            f.write(line)
            f2.write(line)
        else:
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()

            liftpos = lo.convert_coordinate('chr'+arr[0], int(arr[1]))

            ref_b38 = ''
            chrom_b38 = ''
            pos_b38 = ''

            for l1 in range(len(arr[2].strip())):
                for p2 in liftpos:
                    recs = tb_b38.query(p2[0].replace('chr',''), int(p2[1])+l1, int(p2[1])+l1)
                    for r1 in recs:
                        chrom_b38 = p2[0].replace('chr','')
                        pos_b38 = p2[1]
                        ref_b38 += r1[2]
                        break
            if chrom_b38 != '' and arr[2].strip() == ref_b38:
                # print (arr, liftpos, ref_b38)
                # arr[7] += ';GRCh37=' + arr[0] + ':' + arr[1]
                arr[0] = chrom_b38
                arr[1] = str(pos_b38)
                f.write('\t'.join(arr)+'\n')

            else:
                f2.write('\t'.join(arr)+'\n')
            # break
            i += 1
            if i % 10000 == 0:
                print (i, arr)
                # break
                
            pass

    f.close()
    f2.close()

    proc_util.run_cmd('vcf-sort -c ' + out + ' > ' + out + '.sorted.vcf')
    proc_util.run_cmd('tabixgz '+out+ '.sorted.vcf',True)

    proc_util.run_cmd('vcf-sort -c ' + out2 + ' > ' + out2 + '.sorted.vcf')
    proc_util.run_cmd('tabixgz '+out2+ '.sorted.vcf',True)
    

if __name__ == "__main__":
    b38 = "/home/mk446/BiO/Data/CADD/v1.4/GRCh38/whole_genome_SNVs.tsv.gz"
    liftover(sys.argv[1])
