#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### run_insert.py
#### made by Daniel Minseok Kwon
#### 2019-06-05 09:32:32
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

def run_insert():
	#for k in range(2,1000+1):
	for k in range(1001,4000+1):
		print ('./s02_insert.sh ./tmp/gnomad_'+str(k)+'.json')
		# proc_util.run_cmd('./s02_insert.sh ./tmp/gnomad_'+str(k)+'.json', True)
		# proc_util.run_cmd('du -h > rst_insert/du_' + str(k) + '.txt', True)
		# proc_util.run_cmd('elist > rst_insert/elist_' + str(k) + '.txt', True)
		# proc_util.run_cmd('python search_test.py rst_insert/search_test_result_' + str(k) + '.txt', True)

if __name__ == "__main__":
	#c = db_util.connect("")
	run_insert()
