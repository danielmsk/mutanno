#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse
from . import _version

PROG = "mutanno"


def get_options():
    parser = argparse.ArgumentParser(
        usage='%(prog)s <sub-command> [options]', description='%(prog)s ver' + _version.VERSION + " (" +
        _version.VERSION_DATE + ")" + ': python tool for variant annotation')
    parser.add_argument('-v', '--version', action='version',
                        version="%(prog)s ver" + _version.VERSION + " (" + _version.VERSION_DATE + ")")
    subparsers = parser.add_subparsers(
        title="sub-commands", dest="subcommand", metavar='', prog=PROG)

    p1 = subparsers.add_parser('annot', help='annotation', description='annotation')
    p1.add_argument('-vcf', dest='vcf', default='', help='VCF file')
    p1.add_argument('-out', dest='out', default='', help='title of output file')
    p1.add_argument('-ds', dest='ds', default='datastructure.json', help='data structure json file')
    p1.add_argument('-sourcefile', dest='sourcefile', default='', help='data source file')
    p1.add_argument('-remove_unannotated_variant', dest='remove_unannotated_variant', default=False,
                    action="store_true", help='remove unannotated variants in out vcf')
    # p1.add_argument('-buff', dest='buff', type=int, default=100, help='loading size in memory')
    p1.add_argument('-blocksize', dest='blocksize', type=int, default=1000, help='loading size in memory')
    p1.add_argument('-add_genoinfo', dest='add_genoinfo',
                    action="store_true", default=False, help='add genotype info. in INFO field')
    p1.add_argument('-split_multi_allelic_variant', dest='split_multi_allelic_variant',
                    action="store_true", default=False, help='split multi-allelic variants')
    p1.add_argument('-clean_tag', dest='clean_tag_list',
                    default=[], help='remove previous annotation information', nargs='*')
    p1.add_argument('-load_source_in_memory', dest='load_source_in_memory',
                    action="store_true", default=False, help='loading data source in memory')
    p1.add_argument('-sparse', dest='sparse',
                    action="store_true", default=False, help='for sparse variant position')
    p1.add_argument('-log', dest='logfile', default='', help='log file')
    p1.add_argument('-silence', dest='silence', action="store_true",
                    default=False, help='do not print any log.')
    p1.add_argument('-debug', dest='debug', action="store_true",
                    default=False, help='turn on the debugging mode')

    p1 = subparsers.add_parser('makedata', help='make a single data source file',
                               description='make a single data source file')
    p1.add_argument('-out', dest='out', default='', help='title of output file')
    p1.add_argument('-ds', dest='ds', default='', help='datasource json file')
    p1.add_argument('-region', dest='region', default='',
                    help='target region: (ex -region chr1:12345678-22345678 )')
    p1.add_argument('-vartype', dest='vartype', default='all', help='variant type')
    p1.add_argument('-blocksize', dest='blocksize', type=int, default=10000, help='blocksize')
    p1.add_argument('-debug', dest='debug', action="store_true",
                    default=False, help='turn on the debugging mode')
    p1.add_argument('-silence', dest='silence', action="store_true", default=False, help='do not print any log.')

    p1 = subparsers.add_parser('convert', help='convert', description='convert')
    p1.add_argument('-vcf2tsv', dest='vcf2tsv', action="store_true", default=False, help='convert vcf to tsv format')
    p1.add_argument('-vep2tab', dest='vep2tab', action="store_true", default=False, help='convert vep to tsv format')
    p1.add_argument('-in', dest='in', default='', help='input file')
    p1.add_argument('-out', dest='out', default='', help='title of output file')
    p1.add_argument('-ds', dest='ds', default='', help='datastructure json file')
    p1.add_argument('-region', dest='region', default='',
                    help='target region: (ex -region chr1:12345678-22345678 )')
    p1.add_argument('-chromsplit', dest='chromsplit', action="store_true",
                    default=False, help='save separate files by chromosomes')
    p1.add_argument('-silence', dest='silence', action="store_true", default=False, help='do not print any log.')
    # p1.add_argument('-debug', dest='debug', action="store_true", default=False, help='turn on the debugging mode')

    # p1 = subparsers.add_parser('make_datasource', help='convert', description='convert')
    # p1.add_argument('-ds', dest='ds', default='datastructure.json', help='datasource json file')
    # p1.add_argument('-variant', dest='variant', default='datastructure.json', help='datasource json file')

    p1 = subparsers.add_parser('precal', help='pre-calculate', description='pre-calculate')
    p1.add_argument('-make_input_vcf', dest='make_input_vcf', action="store_true",
                    default=False, help='make input VCF file')
    p1.add_argument('-merge_vep', dest='merge_vep',
                    action="store_true", default=False, help='check and merge VEP result')
    p1.add_argument('-check_vep_result', dest='check_vep_result',
                    action="store_true", default=False, help='check VEP result')
    p1.add_argument('-run_vep', dest='run_vep', action="store_true", default=False, help='run vep')
    p1.add_argument('-out', dest='out', default='', help='title of output file')
    p1.add_argument('-fasta', dest='fasta', default='', help='reference')
    p1.add_argument('-vep', dest='vep', default='', help='vep path')
    p1.add_argument('-vepcache', dest='vepcache', default='', help='vep cache directory path')
    p1.add_argument('-cache_version', dest='cache_version', default='98', help='vep cache version')
    p1.add_argument('-region', dest='region', default='', help='region (chr1:123456-789012)')
    p1.add_argument('-vcf', dest='vcf', default="", help='VCF file')
    p1.add_argument('-vep_result', dest='vep_result', default='', help='vep result file')

    p1 = subparsers.add_parser('datacheck', help='check data format',
                               description='check data format')
    # p1.add_argument('-bam', dest='bam', default='', help='BAM file')
    # p1.add_argument('-out', dest='out', default='', help='title of output file')
    # p1.add_argument('-silence', dest='silence', action="store_true", default=False, help='don\'t print any log.')
    # p1.add_argument('-debug', dest='debug', action="store_true", default=False, help='turn on the debugging mode')

    p1 = subparsers.add_parser('preprocess', help='quality metrics for VCF',
                               description='quality metrics for VCF')
    p1.add_argument('-infile', dest='infile', default='', help='title of input file')
    p1.add_argument('-out', dest='out', default='', help='title of output file')
    p1.add_argument('-ds', dest='ds', default='', help='datasource json file')
    p1.add_argument('-make_dbnsfp_transcript', dest='make_dbnsfp_transcript', default=False, action="store_true",
                    help='make dbNSFP transcript file')

    if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1][0] != '-'):
        sys.argv.append('-h')
    opt = vars(parser.parse_args())
    return opt
