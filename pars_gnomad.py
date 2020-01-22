#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### pars_gnomad.py
#### made by Min-Seok Kwon
#### 2019-06-04 12:43:08
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
import seq_util
import tabix
import json
from elasticsearch import Elasticsearch
from datetime import datetime
import time
# es = Elasticsearch()

def pars_gnomad_pyelasticsearch():
	gnomad_vcf = GNOMAD_RAW_VCF.replace('##CHROM##','1')
	tb = tabix.open(gnomad_vcf)
	idx='test-index'
	es.indices.create(index=idx, ignore=400)
	i = 0
	for rec in tb.query("1", 1, 1250000):
		i += 1
		d = {}
		d['CHROM'] = rec[0]
		d['POS'] = int(rec[1])
		d['ID'] = rec[2]
		d['REF'] = rec[3]
		d['ALT'] = rec[4]
		k1 = d['CHROM'] + ':' + str(d['POS']) + '_' + d['REF'] + '>' + d['ALT']

		# body = '{"index":{"_id":"'+k1+'"}}'
		# body += json.dumps(d)
		# print (body)
		# es.bulk(index=idx,body=body)
		d["timestamp"] = datetime.now()

		# print (d)
		print (es.index(index=idx, id=k1, body=d))
		print (es.get(index=idx, id=k1)['_source'])
		if i > 100:
			break

	res = es.search(index="test-index", body={"query": {"match_all": {}}})
	print("Got %d Hits:" % res['hits']['total']['value'])


def convert_gnomad_to_dict(rec):
	d = {}
	d['CHROM'] = rec[0]
	d['POS'] = int(rec[1])
	d['ID'] = rec[2]
	d['REF'] = rec[3]
	d['ALT'] = rec[4]
	d['QUAL'] = rec[5]

	d2 = {}
	no_filed = 0
	for f1 in rec[7].strip().split(';'):
		arr = f1.split('=')
		if len(arr) == 2:
			if arr[0].strip() != "":
				d2[arr[0].strip()] = arr[1].strip()
				no_filed+=1
		else:
			if f1.strip() != "":
				d2[f1.strip()] = "true"
				no_filed+=1
		
	d['1000G_genome'] = d2
	return d


def pars_gnomad():
	idx='test'
	# cmd = 'curl -X PUT "localhost:9200/'+idx+'?pretty"'
	# cmd = 'curl -X PUT "localhost:9200/'+idx+'?pretty" -H "Content-Type: application/json" --data-binary @gnomad_mapping.json'
	# print (proc_util.run_cmd(cmd, True))
	# cmd = 'curl -X PUT "localhost:9200/'+idx+'/_settings" -H "Content-Type: application/json" -d \'{"index": {"mapping.total_fields.limit": 1000000} }\''
	# print (proc_util.run_cmd(cmd, True))
	# cmd = 'curl -X GET "localhost:9200/'+idx+'?pretty"'
	# print (proc_util.run_cmd(cmd, True))


	jsonlist = []
	k = 1
	out = 'tmp/gnomad_'+str(k)+'.json'
	f = open(out,'w')
	jsonlist.append(out)
	for chrom in seq_util.MAIN_CHROM_LIST:
		gnomad_vcf = GNOMAD_RAW_VCF.replace('##CHROM##',chrom)
		tb = tabix.open(gnomad_vcf)
		
		i = 0
		for rec in tb.query(chrom, 1, seq_util.CHROM_LEN['b37d5'][chrom]+1):
			i += 1
			d = convert_gnomad_to_dict(rec)
			k1 = d['CHROM'] + ':' + str(d['POS']) + '_' + d['REF'] + '>' + d['ALT']
			d["timestamp"] = str(datetime.now())
			# print (len(d.keys()))
			body = '{"index":{"_id":"'+k1+'"}}\n'
			body += json.dumps(d)+'\n'
			f.write(body)
			if i % 3000 == 0:
				# break
				f.close()
				k += 1
				out = 'tmp/gnomad_'+str(k)+'.json'
				jsonlist.append(out)
				f = open(out,'w')
				print (chrom, d['POS'], i, out)
		# break
	f.close()
	print ("Saved",out)

	# time.sleep(10)


	for jsonfile in jsonlist:
	# curl -X PUT "localhost:9200/test-index?pretty"
	# cmd = 'curl -X POST "localhost:9200/test-index/_bulk?pretty" -H "Content-Type: application/json" --data-binary @gnomad.json'
		cmd = 'curl -X POST "localhost:9200/'+idx+'/_bulk?pretty&refresh" -H "Content-Type: application/json" --data-binary @'+jsonfile
		print (cmd)
		# proc_util.run_cmd(cmd, True)
	# print (proc_util.run_cmd(cmd, True))

if __name__ == "__main__":
	# GNOMAD_RAW_VCF = "/home/mk446/park/DATA/gnomAD/2.1.1/gnomad.exomes.r2.1.1.sites.##CHROM##.vcf.bgz"
	GNOMAD_RAW_VCF = "/home/mk446/park/DATA/gnomAD/2.1.1/gnomad.genomes.r2.1.1.sites.##CHROM##.vcf.bgz"
	pars_gnomad()
