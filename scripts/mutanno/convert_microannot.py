#!/usr/bin/env python
# -*- coding: utf-8 -*-
# convert_microannot.py
# made by Daniel Minseok Kwon
# 2019-12-18 14:37:51
#########################
import os
import sys
import tabix
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path = "/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path = "/ms1/bin/python_lib"
else:
    sys_path = "/home/mk446/bin/python_lib"
sys.path.append(sys_path)


def get_varmap(fp, k1, chrom, pos):
    varmap = {}
    varmap['pos'] = pos
    varmap['recs'] = {}
    recs = fp[k1].query(chrom, pos, pos)
    for r1 in recs:
        if k1 == 'spliceai':
            ref = r1[3]
            alt = r1[4]
        else:
            ref = r1[2]
            alt = r1[3]
        varmap['recs'][ref + '/' + alt] = r1
    return varmap


def convert_microannot():
    fp = {}
    varmap = {}
    for k1 in fn.keys():
        fp[k1] = tabix.open(fn[k1])
        varmap[k1] = {}
        # print (k1, fp[k1].readline())

    out = "/home/mk446/bio/mutanno/DATASOURCE/MICROANNOT/microannot_datasource_v1.0.tsv"
    out = "microannot_datasource_v1.0.tsv"
    f = open(out, 'w')

    for chrom in seq_util.MAIN_CHROM_LIST:
        for line in file_util.gzopen(vep.replace('#CHROM#', chrom)):
            line = line.decode('UTF-8')
            if line[0] != '#':
                arr = line.split('\t')
                pos = int(arr[1])
                ref = arr[3]
                alt = arr[4]
                varfunc = arr[10]
                gene = arr[22]
                if gene == "-":
                    gene = ""
                info = ["VEP=" + varfunc + "|" + gene]
                cont = arr[:2]
                cont.append(ref)
                cont.append(alt)

                for k1 in fn.keys():
                    if 'pos' not in varmap[k1].keys() or varmap[k1]['pos'] != pos:
                        varmap[k1] = get_varmap(fp, k1, chrom, pos)
                    try:
                        annot = varmap[k1]['recs'][ref + '/' + alt]
                        print(annot)
                        info.append(k1 + '=')
                    except KeyError:
                        pass

                cont.append(';'.join(info))
                print(cont)
                f.write('\t'.join(cont) + '\n')
        break
    f.close()


if __name__ == "__main__":
    import file_util
    import seq_util

    vep = "/home/mk446/bio/mutanno/DATASOURCE/ANNOT/VEP/hg38/b38_WGSNV.vep.chr#CHROM#.vcf.gz"
    dspath = "/home/mk446/bio/mutanno/DATASOURCE/"
    fn = {}
    fn['SpliceAI'] = dspath + "SPLICING/SpliceAI/hg38/spliceai_scores.raw.snv.hg38.vcf.gz"
    fn['gnomADgnome'] = dspath + "POPULATION_AF/gnomAD/hg38/gnomad.genomes.r3.0.sites.AF_popmax.tsv.gz"
    fn['gnomADexome'] = dspath + \
        "POPULATION_AF/gnomAD/hg38/r2.1.1/gnomad.exomes.r2.1.1.sites.liftover_grch38.AF_popmax.tsv.gz"
    fn['CLINVAR'] = dspath + "VARIANTDB/CLINVAR/hg38/clinvar_20191202.vcf.gz"
    convert_microannot()
