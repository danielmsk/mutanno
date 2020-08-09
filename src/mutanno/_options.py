#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse
from . import _version
from .util import file_util

PROG = "mutanno"

DEFAULT_OPT = {'annot': {}, 'makedata': {}, 'convert': {}}
DEFAULT_OPT['annot']['vcf'] = ""
DEFAULT_OPT['annot']['out'] = "mutanno"
DEFAULT_OPT['annot']['outtype'] = ['vcf']
DEFAULT_OPT['annot']['ds'] = ""
DEFAULT_OPT['annot']['sourcefile'] = ""
DEFAULT_OPT['annot']['genoinfo'] = False
DEFAULT_OPT['annot']['hgvs'] = False
DEFAULT_OPT['annot']['variant_class'] = False
DEFAULT_OPT['annot']['hg19'] = False
DEFAULT_OPT['annot']['genetable'] = False
DEFAULT_OPT['annot']['split_multi_allelic_variant'] = False
DEFAULT_OPT['annot']['single_source_mode'] = False
DEFAULT_OPT['annot']['load_source_in_memory'] = False
DEFAULT_OPT['annot']['sparse'] = False
DEFAULT_OPT['annot']['silence'] = False
DEFAULT_OPT['annot']['debug'] = False


def convert_valuetype(typestr):
    rsttype = None
    if typestr is not None:
        if typestr == "int":
            rsttype = int
        if typestr == "float":
            rsttype = float
    return rsttype


def _get_options():
    global DEFAULT_OPT

    confjson = file_util.load_json(file_util.getDataPath('conf.json'))

    parser = argparse.ArgumentParser(usage='%(prog)s <sub-command> [options]',
                                     description='%(prog)s ver' + confjson['VERSION'] +
                                     " (" + confjson['VERSION_DATE'] + ")" + ': convert bam to image')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s ver' + confjson['VERSION'] + " (" + confjson['VERSION_DATE'] + ")")
    parser.add_argument('-silence', dest='silence', action="store_true", default=False, help='don\'t print any log.')

    subparsers = parser.add_subparsers(title="sub-commands", dest="subcommand", metavar='', prog=PROG)

    for s1 in confjson['options']:
        p1 = subparsers.add_parser(s1['subcommand'], help=s1['help'], description=s1['description'])
        for a1 in s1['options']:
            valuetype = convert_valuetype(a1['type'])
            if "dest" in a1.keys() and a1['dest'] != "":
                dest = a1['dest']
            else:
                dest = a1['param']
            if a1['action'] is not None:
                p1.add_argument('-' + a1['param'], dest=dest, default=a1['default'],
                                help=a1['help'], action=a1['action'])
            else:
                p1.add_argument('-' + a1['param'], dest=dest, default=a1['default'],
                                help=a1['help'], nargs=a1['nargs'], type=valuetype)

    '''
    do = DEFAULT_OPT['makedata']
    p1 = subparsers.add_parser('makedata', help='make a single data source file',
                               description='make a single data source file')
    p1.add_argument('-out', dest='out', default='', help='title of output file')
    p1.add_argument('-outtype', dest='outtype', default='json',
                    choices=['tsv', 'json'], help='output file type (default: tsv)')
    p1.add_argument('-ds', dest='ds', default='', help='datasource json file')
    p1.add_argument('-region', dest='region', default='',
                    help='target region: (ex -region chr1:12345678-22345678 )')
    p1.add_argument('-region_vcf', dest='region_vcf', default='',
                    help='target regions using vcf file')
    p1.add_argument('-target_source', dest='target_source', default='', help='target source')
    p1.add_argument('-blocksize', dest='blocksize', default='', help='block size')
    p1.add_argument('-vartype', dest='vartype', default='all',
                    choices=['SNV', 'GENE', 'GENE_MAIN_CHROM', 'CODING_GENE', 'CODING_GENE_MAIN_CHROM'],
                    help='variant type')
    p1.add_argument('-apply_datastructure', dest='apply_datastructure', action="store_true",
                    default=False, help='turn on apply_datastructure mode')
    p1.add_argument('-debug', dest='debug', action="store_true",
                    default=False, help='turn on the debugging mode')

    p1.add_argument('-check', action="store_true", dest='check', default=False, help='check output file')
    p1.add_argument('-log', dest='logfile', default='', help='log file')
    p1.add_argument('-silence', dest='silence', action="store_true", default=False, help='do not print any log.')

    do = DEFAULT_OPT['convert']
    p1 = subparsers.add_parser('convert', help='convert', description='convert')
    p1.add_argument('-vcf2tsv', dest='vcf2tsv', action="store_true", default=False, help='convert vcf to tsv format')
    p1.add_argument('-vcf2html', dest='vcf2html', action="store_true", default=False, help='convert vcf to html format')
    p1.add_argument('-vcf2json', dest='vcf2json', action="store_true", default=False, help='convert vcf to json format')
    p1.add_argument('-tsi2json', dest='tsi2json', action="store_true", default=False, help='convert vcf to json format')
    p1.add_argument('-in', dest='infile', default='', help='input file')
    p1.add_argument('-out', dest='outfile', default='', help='title of output file')
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

    p1 = subparsers.add_parser('validate', help='validate data format',
                               description='validate data format')
    p1.add_argument('-vcf', dest='vcf', default='', help='VCF file')
    p1.add_argument('-ds', dest='ds', default="", help='data structure json file')
    # p1.add_argument('-out', dest='out', default='', help='title of output file')
    p1.add_argument('-silence', dest='silence', action="store_true", default=False, help='don\'t print any log.')
    p1.add_argument('-debug', dest='debug', action="store_true", default=False, help='turn on the debugging mode')

    p1 = subparsers.add_parser('preprocess', help='quality metrics for VCF',
                               description='quality metrics for VCF')
    p1.add_argument('-infile', dest='infile', default='', help='title of input file')
    p1.add_argument('-out', dest='out', default='', help='title of output file')
    p1.add_argument('-ds', dest='ds', default='', help='datasource json file')
    p1.add_argument('-make_dbnsfp_transcript', dest='make_dbnsfp_transcript', default=False, action="store_true",
                    help='make dbNSFP transcript file')

    p1 = subparsers.add_parser('web', help='web view',
                               description='web view')
    p1.add_argument('-ds', dest='ds', default='', help='datasource json file')
    p1.add_argument('-port', dest='port', default='8080', help='port')
    p1.add_argument('-log', dest='logfile', default='', help='log file')
    p1.add_argument('-silence', dest='silence', action="store_true", default=False, help='do not print any log.')
    p1.add_argument('-debug', dest='debug', action="store_true", default=False, help='turn on the debugging mode')

    p1 = subparsers.add_parser('view', help='view', description='view for annotated VCF and tsi file')
    p1.add_argument('-ds', dest='ds', default='', help='datasource json file')
    p1.add_argument('-vcf', dest='vcf', default='', help='VCF file')
    p1.add_argument('-region', dest='region', default='',
                    help='target region: (ex -region chr1:12345678-22345678 )')
    p1.add_argument('-log', dest='logfile', default='', help='log file')
    p1.add_argument('-silence', dest='silence', action="store_true", default=False, help='do not print any log.')
    p1.add_argument('-debug', dest='debug', action="store_true", default=False, help='turn on the debugging mode')
    '''

    if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1][0] != '-'):
        sys.argv.append('-h')
    opt = parser.parse_args()
    # opt = vars(parser.parse_args())
    return opt


def get_options():
    global DEFAULT_OPT

    parser = argparse.ArgumentParser(
        usage='%(prog)s <sub-command> [options]', description='%(prog)s ver' + _version.VERSION + " (" +
        _version.VERSION_DATE + ")" + ': python tool for variant annotation')
    parser.add_argument('-v', '--version', action='version',
                        version="%(prog)s ver" + _version.VERSION + " (" + _version.VERSION_DATE + ")")
    subparsers = parser.add_subparsers(
        title="sub-commands", dest="subcommand", metavar='', prog=PROG)

    do = DEFAULT_OPT['annot']
    p1 = subparsers.add_parser('annot', help='annotation', description='annotation')
    p1.add_argument('-vcf', dest='vcf', default=do['vcf'], help='VCF file')
    p1.add_argument('-out', dest='out', default=do['out'], help='title of output file')
    p1.add_argument('-outtype', dest='outtype', default=do['outtype'],
                    choices=['vcf', 'json'], help='output file type', nargs='*')
    p1.add_argument('-ds', dest='ds', default=do['ds'], help='data structure json file')
    p1.add_argument('-sourcefile', dest='sourcefile', default=do['sourcefile'], help='data source file')
    p1.add_argument('-genoinfo', dest='genoinfo',
                    default=do['genoinfo'], help='add genotype info. in INFO field', nargs="*")
    p1.add_argument('-hgvs', dest='hgvs',
                    action="store_true", default=do['hgvs'], help='add hgvs')
    p1.add_argument('-variant_class', dest='variant_class',
                    action="store_true", default=do['variant_class'], help='add variant class')
    p1.add_argument('-hg19', dest='hg19',
                    action="store_true", default=do['hg19'], help='add hg19 coordinates')
    p1.add_argument('-chain', dest='chain',
                    default="", help='chain file for liftover of hg19 coordinates')
    p1.add_argument('-genetable', dest='add_genetable',
                    action="store_true", default=do['genetable'], help='add gene table')
    p1.add_argument('-blocksize', dest='blocksize', default='', help='block size')
    p1.add_argument('-split_multi_allelic_variant', dest='split_multi_allelic_variant',
                    action="store_true", default=do['split_multi_allelic_variant'], help='split multi-allelic variants')
    p1.add_argument('-fixpl', dest='fixpl',
                    action="store_true", default=False, help='fix PL values for split multi-allelic variants')
    p1.add_argument('-clean_tag', dest='clean_tag_list',
                    default=[], help='remove previous annotation information', nargs='*')
    p1.add_argument('-single_source_mode', dest='single_source_mode',
                    action="store_true", default=do['single_source_mode'], help='single source mode')
    p1.add_argument('-load_source_in_memory', dest='load_source_in_memory',
                    action="store_true", default=do['load_source_in_memory'], help='loading data source in memory')
    p1.add_argument('-sparse', dest='sparse',
                    action="store_true", default=do['sparse'], help='for sparse variant position')
    p1.add_argument('-log', dest='logfile', default='', help='log file')
    p1.add_argument('-silence', dest='silence', action="store_true",
                    default=do['silence'], help='do not print any log.')
    p1.add_argument('-debug', dest='debug', action="store_true",
                    default=do['debug'], help='turn on the debugging mode')

    p1 = subparsers.add_parser('download', help='download data source', description='download data source')
    p1.add_argument('-source_path', dest='source_path', default="", help='source path')
    p1.add_argument('-source', dest="", default="", help='source name')
    p1.add_argument('-version', dest="lastest", default="", help='source version')
    p1.add_argument('-refversion', dest="hg38", default="", help='reference version')

    do = DEFAULT_OPT['makedata']
    p1 = subparsers.add_parser('makedata', help='make a single data source file',
                               description='make a single data source file')
    p1.add_argument('-out', dest='out', default='', help='title of output file')
    p1.add_argument('-outtype', dest='outtype', default='json',
                    choices=['tsv', 'json'], help='output file type (default: tsv)')
    p1.add_argument('-ds', dest='ds', default='', help='datasource json file')
    p1.add_argument('-region', dest='region', default='',
                    help='target region: (ex -region chr1:12345678-22345678 )')
    p1.add_argument('-region_vcf', dest='region_vcf', default='',
                    help='target regions using vcf file')
    p1.add_argument('-target_source', dest='target_source', default='', help='target source')
    p1.add_argument('-blocksize', dest='blocksize', default='', help='block size')
    p1.add_argument('-vartype', dest='vartype', default='all',
                    choices=['SNV', 'GENE', 'GENE_MAIN_CHROM', 'CODING_GENE', 'CODING_GENE_MAIN_CHROM'],
                    help='variant type')
    p1.add_argument('-apply_datastructure', dest='apply_datastructure', action="store_true",
                    default=False, help='turn on apply_datastructure mode')
    p1.add_argument('-debug', dest='debug', action="store_true",
                    default=False, help='turn on the debugging mode')

    p1.add_argument('-check', action="store_true", dest='check', default=False, help='check output file')
    p1.add_argument('-log', dest='logfile', default='', help='log file')
    p1.add_argument('-silence', dest='silence', action="store_true", default=False, help='do not print any log.')

    do = DEFAULT_OPT['convert']
    p1 = subparsers.add_parser('convert', help='convert', description='convert')
    p1.add_argument('-vcf2tsv', dest='vcf2tsv', action="store_true", default=False, help='convert vcf to tsv format')
    p1.add_argument('-vcf2html', dest='vcf2html', action="store_true", default=False, help='convert vcf to html format')
    p1.add_argument('-vcf2json', dest='vcf2json', action="store_true", default=False, help='convert vcf to json format')
    p1.add_argument('-tsi2json', dest='tsi2json', action="store_true", default=False, help='convert vcf to json format')
    # p1.add_argument('-vep2tab', dest='vep2tab', action="store_true", default=False, help='convert vep to tsv format')
    p1.add_argument('-in', dest='infile', default='', help='input file')
    p1.add_argument('-out', dest='outfile', default='', help='title of output file')
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

    p1 = subparsers.add_parser('validate', help='validate data format',
                               description='validate data format')
    p1.add_argument('-vcf', dest='vcf', default='', help='VCF file')
    p1.add_argument('-ds', dest='ds', default="", help='data structure json file')
    # p1.add_argument('-out', dest='out', default='', help='title of output file')
    p1.add_argument('-silence', dest='silence', action="store_true", default=False, help='don\'t print any log.')
    p1.add_argument('-debug', dest='debug', action="store_true", default=False, help='turn on the debugging mode')

    p1 = subparsers.add_parser('preprocess', help='quality metrics for VCF',
                               description='quality metrics for VCF')
    p1.add_argument('-infile', dest='infile', default='', help='title of input file')
    p1.add_argument('-out', dest='out', default='', help='title of output file')
    p1.add_argument('-ds', dest='ds', default='', help='datasource json file')
    p1.add_argument('-make_dbnsfp_transcript', dest='make_dbnsfp_transcript', default=False, action="store_true",
                    help='make dbNSFP transcript file')

    p1 = subparsers.add_parser('web', help='web view',
                               description='web view')
    p1.add_argument('-ds', dest='ds', default='', help='datasource json file')
    p1.add_argument('-port', dest='port', default='8080', help='port')
    p1.add_argument('-log', dest='logfile', default='', help='log file')
    p1.add_argument('-silence', dest='silence', action="store_true", default=False, help='do not print any log.')
    p1.add_argument('-debug', dest='debug', action="store_true", default=False, help='turn on the debugging mode')

    p1 = subparsers.add_parser('view', help='view', description='view for annotated VCF and tsi file')
    p1.add_argument('-ds', dest='ds', default='', help='datasource json file')
    p1.add_argument('-vcf', dest='vcf', default='', help='VCF file')
    p1.add_argument('-region', dest='region', default='',
                    help='target region: (ex -region chr1:12345678-22345678 )')
    p1.add_argument('-log', dest='logfile', default='', help='log file')
    p1.add_argument('-silence', dest='silence', action="store_true", default=False, help='do not print any log.')
    p1.add_argument('-debug', dest='debug', action="store_true", default=False, help='turn on the debugging mode')

    if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1][0] != '-'):
        sys.argv.append('-h')
    opt = parser.parse_args()
    # opt = vars(parser.parse_args())
    return opt


def get_opt_object_from_dict(optdict={}):
    opt = Option()
    opt.convert_from_dict(optdict)
    return opt


class Option:
    def __init__(self):
        pass

    def convert_from_dict(self, optdict={}):
        for k1, v1 in optdict.items():
            self.__dict__[k1] = v1
