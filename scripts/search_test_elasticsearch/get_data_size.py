#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### get_data_size.py
#### made by Daniel Minseok Kwon
#### 2019-06-05 10:37:10
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

def get_size(rstfile):
	size = ""
	for line in open(rstfile):
		arr = line.split('\t')
		arr[-1] = arr[-1].strip()
		if arr[1] == "./elasticsearch-7.1.1/data":
			size = arr[0].replace('> ','')
			break
	return size

def get_size_from_elist(rstfile):
	d = {}
	for line in open(rstfile):
		if line[:6] != "health":
			line = line.strip().replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ')
			arr = line.replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').split(' ')
			d['doc_cnt'] = arr[6]
			d['del_cnt'] = arr[7]
			d['size'] = arr[8]
			break
	return d

def get_data_size():
	f = open('get_data_size.txt','w')
	for k in range(1,1000+1):
		rstfile = "./rst_insert/du_" + str(k) + ".txt"
		size = get_size(rstfile)
		rstfile = "./rst_insert/elist_" + str(k) + ".txt"
		docs = get_size_from_elist(rstfile)
		print (k, size, docs['doc_cnt'],docs['del_cnt'],docs['size'])
		cont = []
		cont.append(str(k))
		cont.append(docs['doc_cnt'])
		cont.append(size)
		cont.append(docs['size'])
		cont.append(docs['del_cnt'])
		f.write('\t'.join(cont)+'\n')
	f.close()
if __name__ == "__main__":
	#c = db_util.connect("")
	get_data_size()
