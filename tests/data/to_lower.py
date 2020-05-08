#!/usr/bin/env python
# -*- coding: utf-8 -*-
# to_lower.py
# made by Daniel Minseok Kwon
# 2020-04-27 15:11:56
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


def to_lower(infile, outfile):
	f = open(outfile, 'w')
	ds = file_util.jsonOpen(infile)
	
	for si in range(len(ds['gene_source'])):
		for fi in range(len(ds['gene_source'][si]['fields'])):
			f1 = ds['gene_source'][si]['fields'][fi]
			if 'name2' in f1.keys():
				ds['gene_source'][si]['fields'][fi]['name2'] = f1['name2'].lower().replace('-','_')
			else:
				if f1['name'] != f1['name'].lower():
					# f1['name2'] = f1['name'].lower()
					ds['gene_source'][si]['fields'][fi]['name2'] = f1['name'].lower().replace('-','_')
					print(ds['gene_source'][si]['fields'][fi]['name2'])
				

	# print(ds['gene_source'][0]['name'])
	file_util.jsonSave(outfile, ds, 2)
	
	f.close()

	

if __name__ == "__main__":
	import proc_util
	import file_util
	to_lower('datastructure_v0.4.3ds_mvp.json', 'datastructure_v0.4.4ds_mvp.json')
