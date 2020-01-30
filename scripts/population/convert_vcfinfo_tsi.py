#!/usr/bin/env python
# -*- coding: utf-8 -*-
# convert_vcfinfo_tsi.py
# made by Daniel Minseok Kwon
# 2020-01-28 16:13:41
#########################
import sys
import os
import time
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path = "/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path = "/ms1/bin/python_lib"
else:
    sys_path = "/home/mk446/bin/python_lib"
sys.path.append(sys_path)


def convert_vcfinfo_tsi(vcf):
	out = vcf.replace('.vcf.gz','') + '.tsi'
	f = open(out, 'w')
	infoid_list = []
	infoidx = 0
	no = 0

	# for line in open('/home/mk446/mutanno/DATASOURCE/POPULATION_AF/UK10K/hg38/header'):
	# 	if line[0] == "#":
	# 		if line[:len('##INFO=<ID=')] == "##INFO=<ID=":
	# 			infoid_list.append(line.split(',')[0].replace('##INFO=<ID=','').strip())
	# 		if line[:len('#CHROM')] == "#CHROM":
	# 			arr = line.split('\t')
	# 			arr[-1] = arr[-1].strip()
	# 			for i in range(len(arr)):
	# 				if arr[i] == "INFO":
	# 					infoidx = i
	# 					arr[i] = '\t'.join(infoid_list)
	# 			f.write('\t'.join(arr) + '\n')

	for line in file_util.gzopen(vcf):
		no += 1
		line = line.decode('UTF-8')
		if line[0] == "#":
			if line[:len('##INFO=<ID=')] == "##INFO=<ID=":
				infoid_list.append(line.split(',')[0].replace('##INFO=<ID=','').strip())
			if line[:len('#CHROM')] == "#CHROM":
				arr = line.split('\t')
				arr[-1] = arr[-1].strip()
				for i in range(len(arr)):
					if arr[i] == "INFO":
						infoidx = i
						arr[i] = '\t'.join(infoid_list)
				f.write('\t'.join(arr) + '\n')
		else:
			arr = line.split('\t')
			arr[-1] = arr[-1].strip()
			infovaluemap = {}
			arr[0] = arr[0].replace('chr','')
			for f1 in arr[infoidx].split(';'):
				if '=' in f1:
					arr2 = f1.split('=')
					try:
						infovaluemap[infoid_list.index(arr2[0].strip())] = arr2[1].strip()
					except ValueError:
						pass
				else:
					try:
						infovaluemap[infoid_list.index(f1.strip())] = "True"
					except ValueError:
						pass
			infovalues = []
			for i in range(len(infoid_list)):
				try:
					infovalues.append(infovaluemap[i])
				except KeyError:
					infovalues.append('')
			arr[infoidx] = '\t'.join(infovalues)
			f.write('\t'.join(arr) + '\n')
			if no % 1000000 == 0:
				print(no, arr[:2])
	f.close()

	time.sleep(30)

	proc_util.run_cmd('tabixgz ' + out)
    

if __name__ == "__main__":
    import proc_util
    import file_util
    # vcf = "bravo-dbsnp-all.vcf.gz"
    vcf = sys.argv[1]
    convert_vcfinfo_tsi(vcf)
