#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### mk_gnomad_mapping.py
#### made by Min-Seok Kwon
#### 2019-06-04 17:23:05
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
import json

def mk_gnomad_mapping():
	d = {}
	d['mappings'] = {}
	d['mappings']['properties'] = {}
	d['mappings']['properties']['CHROM'] = {"type": "keyword"}
	d['mappings']['properties']['POS'] = {"type": "integer"}
	d['mappings']['properties']['REF'] = {"type": "keyword"}
	d['mappings']['properties']['ALT'] = {"type": "keyword"}
	d['mappings']['properties']['ID'] = {"type": "keyword"}
	d['mappings']['properties']['1000G_genome'] = {"properties" : {}}


	for line in open('field.info.txt'):
		arr = line.split('\t')
		arr[-1] = arr[-1].strip()
		# print (arr)
		ftype = FTYPEMAP[arr[2]]
		if arr[0] == "vep":
			ftype = "text"
		d['mappings']['properties']['1000G_genome']["properties"][arr[0]] = {"type": ftype}
	

	file_util.fileSave('gnomad_mapping.json',json.dumps(d), 'w')

if __name__ == "__main__":
	FTYPEMAP = {'Integer':'integer','Float':'half_float','Flag':'boolean','String':'keyword'}
	mk_gnomad_mapping()
