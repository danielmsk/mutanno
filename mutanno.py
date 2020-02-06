#!/usr/bin/env python
# -*- coding: utf-8 -*-
# mutanno.py
# made by Daniel Minseok Kwon (minseok_kwon@hms.harvard.edu)
# 2019-10-29 10:59:25
#########################
import sys
import argparse
import annot
import makedata
import makegenedata
import convert
import preprocess
import precal

VERSION = "0.2.8"
VERSION_DATE = "2020.01.21"
PROG = "mutanno"

# 0.2.3 : encode 'space' to '%20' (remove blank space)
# 0.2.4 : change type of `dbNSFP SiPhy_29way_pi` to list
# 0.2.5 : merge some dbNSFP fields into transcript table (dbNSFPTranscript)
# 0.2.6 : add makedata
# 0.2.7 : update annot module
# 0.2.8 : 

def get_options():
    parser = argparse.ArgumentParser(
        usage='%(prog)s <sub-command> [options]', description='%(prog)s ver' + VERSION + " (" + VERSION_DATE + ")" + ': python tool for variant annotation')
    parser.add_argument('-v', '--version', action='version',
                        version="%(prog)s ver" + VERSION + " (" + VERSION_DATE + ")")
    subparsers = parser.add_subparsers(
        title="sub-commands", dest="subcommand", metavar='', prog=PROG)

    p1 = subparsers.add_parser('annot', help='annotation', description='annotation')
    p1.add_argument('-vcf', dest='vcf', default='', help='VCF file')
    p1.add_argument('-out', dest='out', default='', help='title of output file')
    p1.add_argument('-ds', dest='ds', default='datastructure.json', help='datasource json file')
    p1.add_argument('-temp', dest='temp', default='qcboard_bamqc.html', help='template html file')
    # p1.add_argument('-buff', dest='buff', type=int, default=100, help='loading size in memory')
    p1.add_argument('-blocksize', dest='blocksize', type=int, default=1000, help='loading size in memory')
    p1.add_argument('-load_source_in_memory', dest='load_source_in_memory',
                    action="store_true", default=False, help='loading data source in memory')
    p1.add_argument('-sparse', dest='sparse',
                    action="store_true", default=False, help='for sparse variant position')
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
    # p1.add_argument('-silence', dest='silence', action="store_true", default=False, help='don\'t print any log.')
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


def cli():
    opt = get_options()
    # print (opt)
    if opt['subcommand'] == 'annot' and 'vcf' in opt.keys() and opt['vcf'] != "":
        annotvcf = annot.AnnotVCF(opt)
        if opt['load_source_in_memory']:
            annotvcf.run_by_loading_source()
        else:
            annotvcf.run()
    if opt['subcommand'] == 'makedata' and 'out' in opt.keys() and opt['out'] != "":
        if opt['vartype'] == 'GENE':
            md = makegenedata.GeneDataSourceFile(opt)
            md.make_single_source_file()
        else:
            md = makedata.DataSourceFile(opt)
            md.make_single_source_file()
    if opt['subcommand'] == 'convert' and 'in' in opt.keys() and opt['in'] != "":
        if 'ds' in opt.keys() and opt['ds'] != "":
            cv = convert.CONVFILE(opt)
            cv.convert()
        if 'vcf2tsv' in opt.keys() and opt['vcf2tsv']:
            cv = convert.VCF2TSV(opt)
            cv.convert()
        if 'vep2tab' in opt.keys() and opt['vep2tab']:
            cv = convert.VEP2TAB(opt)
            cv.convert()
    if opt['subcommand'] == 'liftover':
        av = liftover.LiftOverVCF(opt)
        av.run()
    if opt['subcommand'] == 'precal':
        mp = precal.PreCalculate(opt)
        mp.run()
    if opt['subcommand'] == 'preprocess':
        if opt['make_dbnsfp_transcript']:
            obj = preprocess.MakeDbnsfpTranscript(opt['infile'], opt['out'], opt['ds'], opt)
            obj.run()


if __name__ == "__main__":
    cli()
    # print ("#Usage: mutanno.py [vcf file or annot file] [input file type (vcf/annot, default:vcf)] [seqver(b37(default)/b38)]")

    # f_snpeff = 'DATASOURCE/ANNOT/b38/b38_WGSNV_12.vcf.gz.eff.vcf.gz'
    # tb_snpeff = tabix.open(f_snpeff)
    # tb_annova = tabix.open(f_annova)
    # tb_vep = tabix.open(f_vep)
