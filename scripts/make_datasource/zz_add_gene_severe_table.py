#!/usr/bin/env python
# -*- coding: utf-8 -*-
# zz_add_gene_severe_table.py
# made by Daniel Minseok Kwon
# 2020-05-05 16:39:59
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


VEP_CONSEQUENCE_ORDER = {}
VEP_CONSEQUENCE_ORDER['transcript_ablation'] = 1
VEP_CONSEQUENCE_ORDER['splice_acceptor_variant'] = 2
VEP_CONSEQUENCE_ORDER['splice_donor_variant'] = 3
VEP_CONSEQUENCE_ORDER['stop_gained'] = 4
VEP_CONSEQUENCE_ORDER['frameshift_variant'] = 5
VEP_CONSEQUENCE_ORDER['stop_lost'] = 6
VEP_CONSEQUENCE_ORDER['start_lost'] = 7
VEP_CONSEQUENCE_ORDER['transcript_amplification'] = 8
VEP_CONSEQUENCE_ORDER['inframe_insertion'] = 9
VEP_CONSEQUENCE_ORDER['inframe_deletion'] = 10
VEP_CONSEQUENCE_ORDER['missense_variant'] = 11
VEP_CONSEQUENCE_ORDER['protein_altering_variant'] = 12
VEP_CONSEQUENCE_ORDER['splice_region_variant'] = 13
VEP_CONSEQUENCE_ORDER['incomplete_terminal_codon_variant'] = 14
VEP_CONSEQUENCE_ORDER['start_retained_variant'] = 15
VEP_CONSEQUENCE_ORDER['stop_retained_variant'] = 16
VEP_CONSEQUENCE_ORDER['synonymous_variant'] = 17
VEP_CONSEQUENCE_ORDER['coding_sequence_variant'] = 18
VEP_CONSEQUENCE_ORDER['mature_miRNA_variant'] = 19
VEP_CONSEQUENCE_ORDER['5_prime_UTR_variant'] = 20
VEP_CONSEQUENCE_ORDER['3_prime_UTR_variant'] = 21
VEP_CONSEQUENCE_ORDER['non_coding_transcript_exon_variant'] = 22
VEP_CONSEQUENCE_ORDER['intron_variant'] = 23
VEP_CONSEQUENCE_ORDER['NMD_transcript_variant'] = 24
VEP_CONSEQUENCE_ORDER['non_coding_transcript_variant'] = 25
VEP_CONSEQUENCE_ORDER['upstream_gene_variant'] = 26
VEP_CONSEQUENCE_ORDER['downstream_gene_variant'] = 27
VEP_CONSEQUENCE_ORDER['TFBS_ablation'] = 28
VEP_CONSEQUENCE_ORDER['TFBS_amplification'] = 29
VEP_CONSEQUENCE_ORDER['TF_binding_site_variant'] = 30
VEP_CONSEQUENCE_ORDER['regulatory_region_ablation'] = 31
VEP_CONSEQUENCE_ORDER['regulatory_region_amplification'] = 32
VEP_CONSEQUENCE_ORDER['feature_elongation'] = 33
VEP_CONSEQUENCE_ORDER['regulatory_region_variant'] = 34
VEP_CONSEQUENCE_ORDER['feature_truncation'] = 35
VEP_CONSEQUENCE_ORDER['intergenic_variant'] = 36


def get_corder(cstr):
    corder = 999
    for c1 in cstr.split('~'):
        if corder > VEP_CONSEQUENCE_ORDER[c1]:
            corder = VEP_CONSEQUENCE_ORDER[c1]
    return corder


def zz_add_gene_severe_table(dsfile):

    out = dsfile + '_mod.tsi'
    f = open(out, 'w')

    for line in file_util.gzopen(dsfile):
        line = file_util.decodeb(line)
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        if line[0] == '#':
            vepsec = []
            for s1 in arr[-1].split(';'):
                if 'VEP=' in s1:
                    header = s1.split('=')[1].split('|')
                    s1 += ";GENES=ENSG|MOST_SEVERE_TRANSCRIPT|MOST_SEVERE_CONSEQUENCE"
                vepsec.append(s1)
            arr[-1] = ';'.join(vepsec)
            
        else:
            veplist = []
            for s1 in arr[-1].split(';'):
                if 'VEP=' in s1:
                    for f1 in s1.split('=')[1].split(','):
                        d = {}
                        for j, value in enumerate(f1.split('|')):
                            d[header[j]] = value
                        veplist.append(d)
            # print(len(veplist))

            sever_transcript = {}
            for j, v1 in enumerate(veplist):
                corder = get_corder(v1['Consequence'])
                enst = v1['Feature']
                ensg = v1['Gene']
                v1['corder'] = corder
                # print(corder, v1['Consequence'], v1['MOST_SEVERE'], v1['Gene'], v1['Feature'])
                try:
                    sever_transcript[ensg]
                    if sever_transcript[ensg]['corder'] > v1['corder']:
                        sever_transcript[ensg] = v1
                    elif sever_transcript[ensg]['corder'] == v1['corder']:
                        if v1['CANONICAL'] == "1":
                            sever_transcript[ensg] = v1
                except KeyError:
                    sever_transcript[ensg] = v1

            for j, v1 in enumerate(veplist):
                enst = v1['Feature']
                ensg = v1['Gene']
                if sever_transcript[ensg]['Feature'] == enst:
                    veplist[j]['MOST_SEVERE'] = "1"
                else:
                    veplist[j]['MOST_SEVERE'] = "0"

            newannot = []
            for s1 in arr[-1].split(';'):
                if 'VEP=' in s1:
                    vepsec = []
                    for j, f1 in enumerate(s1.split('=')[1].split(',')):
                        vepannot = f1.split('|')
                        vepannot[header.index("MOST_SEVERE")] = veplist[j]['MOST_SEVERE']
                        vepsec.append('|'.join(vepannot))
                    s1 = 'VEP=' + ','.join(vepsec)

                    vepsec = []
                    for ensg in sever_transcript.keys():
                        # print(ensg)
                        vepannot = []
                        vepannot.append(ensg)
                        vepannot.append(sever_transcript[ensg]['Feature'])
                        vepannot.append(sever_transcript[ensg]['Consequence'])
                        vepsec.append('|'.join(vepannot))
                    s1 += ";GENES=" + ','.join(vepsec)

                newannot.append(s1)

            arr[-1] = ';'.join(newannot)
            # print(arr)
        line = '\t'.join(arr)
        f.write(line + '\n')

    f.close()
    print('Saved', out)


    

if __name__ == "__main__":
    import proc_util
    import file_util
    dsfile = "/home/mk446/bio/mutanno/DATASOURCE/MAINANNOT/mc_3k.tsi"
    # dsfile = "/home/mk446/bio/mutanno/DATASOURCE/MAINANNOT/tmp/mc_1.tsi"
    zz_add_gene_severe_table(dsfile)
