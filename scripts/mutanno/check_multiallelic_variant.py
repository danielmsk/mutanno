#!/usr/bin/env python
# -*- coding: utf-8 -*-
# scripts/mutanno/check_multiallelic_variant.py
# made by Daniel Minseok Kwon
# 2020-02-25 12:22:58
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

def get_gt(gtvcf):
	arr = gtvcf.split(':')
	gt = arr[0]
	ad = arr[1]
	return gt, ad

def check_multiallelic_variant(rawvcf, annotvcf, out):
	samplelist = []
	gtmap = {}
	f = open(out, 'w')

	for line in file_util.gzopen(annotvcf):
		if annotvcf.endswith('.gz'):
			line = line.decode('UTF-8')
		arr = line.split('\t')
		arr[-1] = arr[-1].strip()
		if line[0] == '#':
			if arr[0] == "#CHROM":
				samplelist = arr[9:]
		else:
			k1 = arr[0] + '_' + arr[1]

			gtlist = [arr[3], arr[4]]
			for j in range(9, len(arr)):
				gt, ad = get_gt(arr[j])
				gtlist.append(gt + ':' +ad)

			try:
				gtmap[k1].append(gtlist)
			except KeyError:
				gtmap[k1] = []
				gtmap[k1].append(gtlist)


	for line in file_util.gzopen(rawvcf):
		if rawvcf.endswith('.gz'):
			line = line.decode('UTF-8')
		arr = line.split('\t')
		arr[-1] = arr[-1].strip()
		if line[0] == '#':
			if arr[0] == "#CHROM":
				cont = arr[:2]
				cont.append(arr[3])
				cont.append(arr[4])
				cont.extend(samplelist)
				cont.append(arr[3])
				cont.append(arr[4])
				cont.extend(samplelist)
				f.write('\t'.join(cont) + '\n')
			
		else:
			k1 = arr[0] + '_' + arr[1]

			for annotgtlist in gtmap[k1]:
				cont = arr[:2]
				cont.append(arr[3])
				cont.append(arr[4])
				for j in range(9, len(arr)):
					gt, ad = get_gt(arr[j])
					cont.append(gt + ':' +ad)				
				cont.extend(annotgtlist)
				f.write('\t'.join(cont) + '\n')

	f.close()
	print('Saved', out)

if __name__ == "__main__":
	import proc_util
	import file_util
	rawvcf = sys.argv[1]
	annotvcf = sys.argv[2]
	out = annotvcf + '.multiallelic_check.txt'
	check_multiallelic_variant(rawvcf, annotvcf, out)
