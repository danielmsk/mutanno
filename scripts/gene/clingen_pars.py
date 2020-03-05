#!/usr/bin/env python
# -*- coding: utf-8 -*-
# clingen_pars.py
# made by Daniel Minseok Kwon
# 2020-03-03 15:42:28
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

def pars_clingen_gene_list(cont):
	genelist = []
	scont = str_util.substr(cont, '<tbody>', '</tbody>').strip()
	for block in scont.split('<tr>'):
		block = block.strip()
		if block != "":
			gblock = str_util.substr(block, "/kb/genes/HGNC:", "</a>")
			arr2 = gblock.split('">')
			hgnc_id = arr2[0].strip()
			gene = arr2[1].strip()
			print(hgnc_id, gene)
			if hgnc_id != "":
				genelist.append([hgnc_id, gene])
			# break
	return genelist

def clingen_pars(url):
	timekey = time_util.getToday()
	html = path + "clingen_list_" + timekey + ".html"
	# cont = web_util.get_url(url)
	# file_util.fileSave(html, cont, 'w')

	cont = file_util.fileOpen(html)
	clingen_gene_list = pars_clingen_gene_list(cont)
	out = html + ".list"
	outcont = "HCNC\tGENE\n"
	for g1 in clingen_gene_list:
		outcont += '\t'.join(g1) + '\n'
	file_util.fileSave(out, outcont, 'w')
	print("Saved", out)



if __name__ == "__main__":
	import proc_util
	import file_util
	import web_util
	import time_util
	import str_util
	url = "https://search.clinicalgenome.org/kb/curations"
	path = "/home/mk446/mutanno/DATASOURCE/GENE/CLINGEN/"
	clingen_pars(url)
