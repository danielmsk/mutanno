#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
def test_check_denovo(vcf):
    i = 0
    for line in file_util.gzopen(vcf):
        if vcf.endswith('.vcf.gz'):
            line = line.decode('UTF-8')
        if line[0] != "#":
            arr = line.split('\t')
            gts = arr[9].split(':')[0]
            gtf = arr[10].split(':')[0]
            gtm = arr[11].split(':')[0]
            # print (gts, gtf, gtm)
            vqsr_filter = arr[6]

            if gtf == '0/0' and gtm == '0/0' and gts != '0/0' and gts != './.':
                i += 1
                print(i, arr[:2], gts, gtf, gtm, vqsr_filter)
        elif line[:2] == "#C":
            arr = line.split('\t')
            print(arr[9:12])
'''
