#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### mk_microannot.py
#### made by Min-Seok Kwon
#### 2019-12-17 12:44:18
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

# poplist = ['AF_popmax','AF_afr','AF_amr','AF_asj','AF_eas','AF_fin','AF_nfe','AF_oth']
# poplist = ['AF_afr','AF_amr','AF_eas','AF_nfe','AF_oth']
poplist = ['AF_afr','AF_amr','AF_eas','AF_nfe']

def check_popmax_list():
    i = 0
    for chrom in seq_util.MAIN_CHROM_LIST:
        vcf = pathgenome.replace('#CHROM#',chrom)
        if file_util.is_exist(vcf):
            print (vcf)
            print ('\t'.join(poplist))
            for line in file_util.gzopen(vcf):
                line = line.decode('UTF-8')
                if line[0] != "#":
                    i += 1
                    arr = line.split('\t')
                    # print(arr[7].split(';'))
                    dat = {'AF_popmax':'0.0'}
                    for p1 in poplist:
                        dat[p1] = ''
                    for f1 in arr[7].split(';'):
                        if "=" in f1 and "AF" in f1:
                            arr2 = f1.split('=')
                            dat[arr2[0].strip()] = arr2[1].strip()
                    # print (dat['AF_popmax'])
                    cont = []
                    aflist = []
                    for p1 in poplist:
                        # print (p1, dat[p1])
                        if dat[p1] == "":
                            dat[p1] = "0.0"
                        aflist.append(float(dat[p1]))
                        cont.append(dat[p1])
                    # print (dat['AF_popmax'])
                    cont.append(str(max(aflist)))
                    cont.append(dat['AF_popmax'])

                    if float(dat['AF_popmax']) == max(aflist):
                        flag = 'O'
                    else:
                        flag = 'X'
                    cont.append(flag)

                    if flag=='X':
                        print ('\t'.join(cont))
                    if i > 1000000:
                        break
        break


def s01_mk_popmaxaf(vcf, out):
    f = open(out, 'w')
    cont = ["#CHROM","POS","ID","REF","ALT","AF_popmax"]
    f.write('\t'.join(cont)+'\n')
    i = 0
    # for chrom in seq_util.MAIN_CHROM_LIST:
    if True:
        # vcf = vcfpath.replace('#CHROM#',chrom)
        if file_util.is_exist(vcf):
            print (vcf)
            print ('Saving',out)
            for line in file_util.gzopen(vcf):
                line = line.decode('UTF-8')
                if line[0] != "#":
                    i += 1
                    arr = line.split('\t')
                    # print(arr[7].split(';'))
                    dat = {'AF_popmax':'0.0'}
                    for p1 in poplist:
                        dat[p1] = '0.0'
                    for f1 in arr[7].split(';'):
                        if "=" in f1 and "AF" in f1:
                            arr2 = f1.split('=')
                            dat[arr2[0].strip()] = arr2[1].strip()
                    
                    aflist = []
                    for p1 in poplist:
                        # print (p1, dat[p1])
                        if dat[p1] == "":
                            dat[p1] = "0.0"
                        aflist.append(float(dat[p1]))
                    
                    popmax = max(aflist)
                    if popmax > 0:
                        arr[2] = ""
                        arr[0] = arr[0].replace('chr','')
                        # if ',' in arr[3] or ',' in arr[4]:
                        #   print (arr[:5], popmax)
                        cont = arr[:5]
                        cont.append(str(popmax))
                        f.write('\t'.join(cont)+'\n')

                    if i % 10000 == 0:
                        print (i, arr[:5])

    f.close()
    file_util.fileSave(out+'.done','','w')


def s02_mk_spliceAI(vcf, out):
    cate = {'<0.2':(0.0,0.2),'0.2~0.4':(0.2,0.4),'0.4~0.6':(0.4,0.6),'0.6~0.8':(0.6,0.8),'0.8~1.0':(0.8,1.0)}
    cntcate = {}

    for c1 in cate.keys():
        cntcate[c1] = 0

    f = open(out, 'w')
    i = 0
    cont = ["#CHROM","POS","ID","REF","ALT","DS_AG","DS_AL","DS_DG","DS_DL"]
    f.write('\t'.join(cont)+'\n')
    cnt = 0
    for line in file_util.gzopen(vcf):
        line = line.decode('UTF-8')
        if line[0] != "#":
            i += 1
            arr = line.split('\t')
            sclist = arr[7].strip().split('|')[2:6]
            flag = False
            fsclist = []
            for sc in sclist:
                fsc = float(sc)
                if fsc >= 0.8:
                    flag = True
                fsclist.append(fsc)
            max_fsc = max(fsclist)

            for c1 in cate.keys():
                if max_fsc >= cate[c1][0] and max_fsc < cate[c1][1]:
                    cntcate[c1] += 1
                    break

            if flag:
                # print (arr)
                cnt += 1
                cont = arr[:5]
                cont.extend(sclist)
                f.write('\t'.join(cont)+'\n')
            if i % 100000 == 0:
                print (i, arr[:5], cnt)
                # break
    f.close()


    cont = ""
    for c1 in cate.keys():
        cont += c1 + '\t' + str(cntcate[c1]) + '\n'
    file_util.fileSave(out+'.stat',cont, 'w')



def s03_mk_cadd(vcf, out):
    cate = {'<5':(0,5),'5~10':(5,10),'10~15':(10,15),'15~20':(15,20),'20~25':(20,25),'25~30':(25,30),'>30':(30,999)}
    cntcate = {}

    for c1 in cate.keys():
        cntcate[c1] = 0

    f = open(out, 'w')
    i = 0
    cont = ["#CHROM","POS","REF","ALT","PHRED"]
    f.write('\t'.join(cont)+'\n')
    cnt = 0
    for line in file_util.gzopen(vcf):
        line = line.decode('UTF-8')
        if line[0] != "#":
            i += 1
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            
            flag = False
            sc = float(arr[5])
            if sc >= 20:
                flag = True

            for c1 in cate.keys():
                if sc >= cate[c1][0] and sc < cate[c1][1]:
                    cntcate[c1] += 1
                    break

            if flag:
                cnt += 1
                cont = arr[:4]
                cont.append(arr[5])
                f.write('\t'.join(cont)+'\n')
            if i % 100000 == 0:
                print (i, arr[:5], cnt)
                # break
    f.close()

    cont = ""
    for c1 in cate.keys():
        cont += c1 + '\t' + str(cntcate[c1]) + '\n'
    file_util.fileSave(out+'.stat',cont, 'w')

def run(pathvcf):
    for s1 in seq_util.get_split_region(chunksize=1000000, seqver='b38d'):
        print (s1)
    # for chrom in seq_util.MAIN_CHROM_LIST:
    #     vcf = pathvcf.replace('#CHROM#',chrom)
    #     out = vcf.replace('.vcf.gz','') + '.AF_popmax.vcf'
    #     if file_util.is_exist(vcf):
    #         print('python /home/mk446/bio/mutanno/DATASOURCE/mk_microannot.py ' + vcf + ' ' + out)


if __name__ == "__main__":
    pathgenome = "/home/mk446/bio/mutanno/DATASOURCE/POPULATION_AF/gnomAD/hg38/r2.1.1/gnomad.genomes.r2.1.1.sites.#CHROM#.liftover_grch38.vcf.gz"
    pathexome = "/home/mk446/bio/mutanno/DATASOURCE/POPULATION_AF/gnomAD/hg38/r2.1.1/gnomad.exomes.r2.1.1.sites.#CHROM#.liftover_grch38.vcf.gz"
    pathgenome_hg38 = "/home/mk446/bio/mutanno/DATASOURCE/POPULATION_AF/gnomAD/hg38/gnomad.genomes.r3.0.sites.chr#CHROM#.vcf.gz"

    # run(pathgenome)
    # run(pathexome)
    # run(pathgenome_hg38)
    
    # check_popmax_list()
    # vcf = sys.argv[1]
    # out = sys.argv[2]
    # s01_mk_popmaxaf(vcf, out)

    out = "/home/mk446/bio/mutanno/DATASOURCE/POPULATION_AF/gnomAD/hg38/r2.1.1/gnomad.genomes.r2.1.1.sites.liftover_grch38.AF_popmax.vcf"
    # s01_mk_popmaxaf(pathgenome, out)
    out = "/home/mk446/bio/mutanno/DATASOURCE/POPULATION_AF/gnomAD/hg38/r2.1.1/gnomad.exomes.r2.1.1.sites.liftover_grch38.AF_popmax.vcf"
    # s01_mk_popmaxaf(pathexome, out)

    
    out = "/home/mk446/bio/mutanno/DATASOURCE/POPULATION_AF/gnomAD/hg38/gnomad.genomes.r3.0.sites.AF_popmax.vcf"
    # s01_mk_popmaxaf(pathgenome_hg38, out)


    s02_mk_spliceAI("/home/mk446/bio/mutanno/DATASOURCE/SPLICING/SpliceAI/hg38/spliceai_scores.raw.snv.hg38.vcf.gz", "/home/mk446/bio/mutanno/DATASOURCE/SPLICING/SpliceAI/hg38/spliceai_scores.raw.snv.hg38.over0.8.vcf")


    CADD = "/home/mk446/bio/mutanno/DATASOURCE/PATHOGENICITY/CADD/hg38/v1.5/whole_genome_SNVs.tsv.gz"
    out = "/home/mk446/bio/mutanno/DATASOURCE/PATHOGENICITY/CADD/hg38/v1.5/whole_genome_SNVs_over20.tsv"
    # s03_mk_cadd(CADD, out)

