#!/usr/bin/env python
# -*- coding: utf-8 -*-
# preproc_HGNC.py
# made by Daniel Minseok Kwon
# 2020-02-24 12:20:11
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


def preproc_HGNC(orifile, outfile):
	f = open(outfile, 'w')
	i = 0
	h = {}
	for line in open(orifile):
		arr = line.split('\t')
		arr[-1] = arr[-1].strip()
		if i == 0:
			for k in range(len(arr)):
				h[arr[k]] = k
			arr[0] = "#" + arr[0]
		else:
			arr[h['hgnc_id']] = arr[h['hgnc_id']].replace('HGNC:', '').strip()

		f.write('\t'.join(arr) + '\n')
		i += 1
	f.close()
	

if __name__ == "__main__":
	import proc_util
	import file_util
	orifile = "/home/mk446/mutanno/DATASOURCE/HGNC/hgnc_complete_set.txt"
	outfile = "/home/mk446/mutanno/DATASOURCE/HGNC/hgnc_complete_set.mod.txt"
	preproc_HGNC(orifile, outfile)
