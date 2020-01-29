#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### s03_liftover_merge.py
#### made by Min-Seok Kwon
#### 2019-11-05 15:47:02
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

def s03_liftover_merge(chrom):
	out = "/home/mk446/bio/mutanno/s03_liftover_map/liftover_b37_b38_chr" + chrom + ".tsv"
	out2 = "/home/mk446/bio/mutanno/s03_liftover_map/liftover_b37_b38_chr" + chrom + "_notmatch.tsv"

	f = open(out, 'w')
	f2 = open(out2, 'w')
	flag = True
	init_i = -1
	for i in range(1,10000000000):
		infile = "/home/mk446/bio/mutanno/s03_liftover_map/" + chrom + "/" + chrom + "_" + str(i) + ".tsv"
		if file_util.is_exist(infile):
			print (infile)
			init_i = i
			for line in open(infile):
				if line[0] == "#" and flag:
					arr = line.split('\t')
					arr[-1] = arr[-1].strip()
					arr.append('MATCH_TYPE')
					f.write('\t'.join(arr)+'\n')
					f2.write('\t'.join(arr)+'\n')
					flag = False
				if line[0] != "#":

					arr = line.split('\t')
					arr[-1] = arr[-1].strip()
					mtype = "."
					if arr[-1] == "":
						mtype = "N"
					if  ";" in arr[-1]:
						mtype = "M"
					if  arr[2] != arr[-1].split('_')[-1]:
						mtype = "G"
					arr.append(mtype)
					if mtype != '.':
						f2.write('\t'.join(arr)+'\n')
					f.write('\t'.join(arr)+'\n')
		else:
			if init_i > 0:
				break
	f.close()
	f2.close()

	proc_util.run_cmd('tabixgz '+out, True)
	proc_util.run_cmd('tabixgz '+out2, True)

def check_files(chrom):
	init_i = -1
	for i in range(1,10000000000):
		infile = "/home/mk446/bio/mutanno/s03_liftover_map/" + chrom + "/" + chrom + "_" + str(i) + ".tsv"
		if file_util.is_exist(infile):
			init_i = i
		else:
			if init_i > 0:
				print (infile)


if __name__ == "__main__":
	chrom = sys.argv[1]
	s03_liftover_merge(chrom)
	# check_files(chrom)
