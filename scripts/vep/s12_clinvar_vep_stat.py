#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s12_clinvar_vep_stat.py
# made by Daniel Minseok Kwon
# 2020-02-10 15:29:47
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


def s12_clinvar_vep_stat():
	

if __name__ == "__main__":
	import proc_util
	import file_util
	s12_clinvar_vep_stat()
