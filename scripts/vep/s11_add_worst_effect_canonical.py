#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s11_add_worst_effect_canonical.py
# made by Daniel Minseok Kwon
# 2020-05-04 09:06:22
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
sys.path.append('..')


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

VEP_BIOTYPE_ORDER = {}
VEP_BIOTYPE_ORDER['protein_coding'] = 1
VEP_BIOTYPE_ORDER['processed_transcript'] = 5
VEP_BIOTYPE_ORDER['transcribed_processed_pseudogene'] = 10
VEP_BIOTYPE_ORDER['processed_pseudogene'] = 12
VEP_BIOTYPE_ORDER['transcribed_unprocessed_pseudogene'] = 15
VEP_BIOTYPE_ORDER['unprocessed_pseudogene'] = 20
VEP_BIOTYPE_ORDER['miRNA'] = 30
VEP_BIOTYPE_ORDER['snRNA'] = 35
VEP_BIOTYPE_ORDER['lncRNA'] = 40
VEP_BIOTYPE_ORDER['enhancer'] = 50
VEP_BIOTYPE_ORDER['promoter'] = 55
VEP_BIOTYPE_ORDER['promoter_flanking_region'] = 60
VEP_BIOTYPE_ORDER['CTCF_binding_site'] = 70
VEP_BIOTYPE_ORDER['TF_binding_site'] = 80
VEP_BIOTYPE_ORDER['open_chromatin_region'] = 99
VEP_BIOTYPE_ORDER[''] = 100

def s11_add_worst_effect_canonical(veptsi):
    
    for chrom in seq_util.MAIN_CHROM_LIST:
        bt = {}
        h = {}

        i = 0
        for line in file_util.gzopen(veptsi.replace('#CHROM#', chrom)):
            line = file_util.decodeb(line)
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()

            if line[0] == "#":
                if '#CHROM' in line:
                    arr[0] = arr[0][1:]
                    arr2 = arr[5].split('|')
                    for k in range(len(arr2)):
                        h[arr2[k]] = k

            else:
                i += 1
                b = []
                g = []

                cano = {}
                
                worst_transcript = {}
                transcripts = []
                gene_level_most_severe = {}
                transcript_level_most_severe = {}
                for f1 in arr[5].split(','):
                    arr2 = f1.split('|')
                    d = {}
                    d['Feature'] = arr2[h['Feature']]
                    d['Gene'] = arr2[h['Gene']]
                    d['SYMBOL'] = arr2[h['SYMBOL']]
                    d['Consequence'] = arr2[h['Consequence']]
                    d['Consequence_order'] = VEP_CONSEQUENCE_ORDER[arr2[h['Consequence']]]
                    d['CANONICAL'] = arr2[h['CANONICAL']]
                    d['BIOTYPE'] = arr2[h['BIOTYPE']]
                    d['BIOTYPE_ORDER'] = VEP_BIOTYPE_ORDER[d['BIOTYPE']]

                    try:
                        transcript_level_most_severe[d['Gene']]
                        if transcript_level_most_severe[d['Gene']]['Consequence_order'] > d['Consequence_order']:
                            transcript_level_most_severe[d['Gene']] = d
                        elif transcript_level_most_severe[d['Gene']]['Consequence_order'] == d['Consequence_order']:
                            if transcript_level_most_severe[d['Gene']]['CANONICAL'] != "YES" and d['CANONICAL'] == "YES":
                                transcript_level_most_severe[d['Gene']] = d
                    except KeyError:
                        transcript_level_most_severe[d['Gene']] = d

                    try:
                        if gene_level_most_severe['Consequence_order'] > d['Consequence_order']:
                            gene_level_most_severe = d
                        elif gene_level_most_severe['Consequence_order'] == d['Consequence_order']:
                            if gene_level_most_severe['BIOTYPE_ORDER'] > d['BIOTYPE_ORDER']:
                                gene_level_most_severe = d
                            elif gene_level_most_severe['BIOTYPE_ORDER'] == d['BIOTYPE_ORDER']:
                                if gene_level_most_severe['CANONICAL'] != "YES" and d['CANONICAL'] == "YES":
                                    gene_level_most_severe = d
                    except KeyError:
                        gene_level_most_severe = d

                    try:
                        bt[d['BIOTYPE']] += 1
                    except KeyError:
                        bt[d['BIOTYPE']] = 1
                    g.append(d['Gene'])
                    b.append(d['BIOTYPE'])
                    transcripts.append(d)
                    print(d)

                print('==>',gene_level_most_severe)
                print('==>',transcript_level_most_severe)
                print()

                # print(transcripts)
                # print(g, b)
                # print()
                # break
                if i % 10 == 0:
                    # for t1 in transcripts:
                    #     print(t1)
                    print()
                    break
        # print(bt)
        break
    

if __name__ == "__main__":
    import file_util
    import seq_util
    import preproc_util
    veptsi = preproc_util.DATASOURCEPATH + '/ANNOT/VEP/hg38/v99/vep.99.hg38.#CHROM#.sorted.rmcsq.tsi.gz'
    s11_add_worst_effect_canonical(veptsi)
