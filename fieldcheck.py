#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### fieldcheck.py
#### made by Min-Seok Kwon
#### 2019-11-21 15:52:07
#########################
import sys
import os

import file_util
import proc_util
import seq_util

def fieldcheck():
	stat = {}
	stat['maxlen'] = {}
	stat['maxex'] = {}
	stat['max_list_size'] = {}
	stat['list_maxlen'] = {}
	stat['list_maxex'] = {}
	for chrom in seq_util.MAIN_CHROM_LIST:
		fname2 = fname.replace('#CHROM#',chrom)
		if file_util.is_exist(fname2):
			print (fname2)
			header = []
			line_no = 1
			for line in file_util.gzopen(fname2):
				line = line.decode('UTF-8')
				if line[0] == '#':
					if line[:len('#C')] == "#C":
						header = line[1:].split('\t')
						header[-1] = header[-1].strip()
				else:
					line_no += 1
					arr = line.split('\t')
					arr[-1] = arr[-1].strip()

					for i in range(len(header)):
						h1 = header[i]
						v1 = arr[i]
						if v1 == '-':
							v1 = ''

						if ',' in arr[i]:
							# print (h1, v1)
							arrv1 = v1.split(',')

							try:
								if stat['max_list_size'][h1] < len(arrv1):
									stat['max_list_size'][h1] = len(arrv1)
							except KeyError:
								stat['max_list_size'][h1] = len(arrv1)

							for v2 in arrv1:
								try:
									if stat['list_maxlen'][h1] < len(v2):
										stat['list_maxlen'][h1] = len(v2)
										stat['list_maxex'][h1] = v2
								except KeyError:
									stat['list_maxlen'][h1] = len(v2)
									stat['list_maxex'][h1] = v2

						try:
							if stat['maxlen'][h1] < len(v1):
								stat['maxlen'][h1] = len(v1)
								stat['maxex'][h1] = v1
						except KeyError:
							stat['maxlen'][h1] = len(v1)
							stat['maxex'][h1] = v1
					# print (line_no, arr[:6])
					if line_no % 10000000 == 0:
						break


			print_stat(stat)
			break


def print_stat(stat):
	for k1 in stat.keys():
		print (k1)
		for k2 in stat[k1].keys():
			cont = "\t"+k2 + " : " + str(stat[k1][k2])
			print (cont)

if __name__ == "__main__":
	print ("#USAGE: python fieldcheck.py [SOURCE] ")
	source = "VEP"
	fname = "/home/mk446/bio/mutanno/DATASOURCE/ANNOT/VEP/hg38/b38_WGSNV.vep.chr#CHROM#.vcf.gz"
	fieldcheck()
