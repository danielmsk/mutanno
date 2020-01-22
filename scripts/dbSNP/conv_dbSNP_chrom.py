#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### conv_dbSNP_chrom.py
#### made by Min-Seok Kwon
#### 2019-11-05 13:50:14
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

def mk_chrom_map(vcf):
	out = ""
	chrommap = {}
	i = 0
	for line in file_util.gzopen(vcf):
		line = line.decode('UTF-8')
		if line[0] == "#":
			pass
		else:
			arr = line.split('\t')
			try:
				chrommap[arr[0]] += 1
			except KeyError:
				chrommap[arr[0]] = 1
			i += 1
			if i % 1000000 == 0:
				print (i, arr[:5])
				# break

	print (chrommap)

	cont = ''
	for chrom in chrommap.keys():
		cont += chrom + '\t' + str(chrommap[chrom]) + '\n'	
	file_util.fileSave("dbSNP_chrom_map.txt", cont, 'w')

def load_chrommap():
	chrommap = {}
	for line in open('dbSNP_chrom_map.txt'):
		arr = line.split('\t')
		if arr[2].strip() != "":
			chrommap[arr[0].strip()] = arr[2].strip()
	return chrommap

def conv_dbSNP_chrom(vcf):
	chrommap = load_chrommap()

	out = vcf.replace('.vcf.gz','.convchrom.vcf')
	i = 0
	f = open(out,'w')
	for line in file_util.gzopen(vcf):
		line = line.decode('UTF-8')
		if line[0] == "#":
			pass
			f.write(line)
		else:
			arr = line.split('\t')
			try:
				arr[0] = chrommap[arr[0]]
				f.write('\t'.join(arr))	
			except KeyError:
				break
			i += 1
			if i % 1000000 == 0:
				print (i, arr[:5])
	f.close()
	proc_util.run_cmd('tabixgz ' + out)

if __name__ == "__main__":
	vcf = "GRCh38_latest_dbSNP_all.vcf.gz"
	# mk_chrom_map(vcf)
	conv_dbSNP_chrom(vcf)
