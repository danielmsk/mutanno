#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s02_check_between_chunk.py
# made by Daniel Minseok Kwon
# 2020-05-05 05:45:09
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


def s02_check_between_chunk(outfile, chunksize):

	regionlist = seq_util.get_split_region(chunksize)

	spos = -9
	for r1 in regionlist:
		k = r1[3]
		out2 = outfile.replace('##',str(k))
		if not file_util.is_exist(out2):
			print('File is not exist.' + out2)
			# cmd = "mv /home/mk446/mutanno/DATASOURCE/MICROANNOT/tmp_splitrun.sh_"+str(k)+".sh /home/mk446/jobs/."
			# proc_util.run_cmd(cmd)
		else:
			cmd = "head -2 " + out2
			cont = proc_util.run_cmd(cmd).strip().split('\n')
			if len(cont) == 2:
				arr = cont[1].split('\t')
				chrom = arr[0]
				spos = int(arr[1])
				
				cmd = "tail -1 " + out2
				cont = proc_util.run_cmd(cmd).strip().split('\n')
				arr = cont[0].split('\t')
				epos = int(arr[1])

				if r1[2] != epos:
					print(epos-r1[2], epos, r1, out2)

			# print(out2)
			

	

if __name__ == "__main__":
	import proc_util
	import file_util
	import seq_util
	# total_chunk_num = sys.argv[1]
	# total_chunk_num = sys.argv[1]
	outfile = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/tmp/mc_##.tsv.tsi"
	chunksize = 1000000
	s02_check_between_chunk(outfile, chunksize)
