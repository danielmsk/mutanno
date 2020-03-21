#!/usr/bin/env python
# -*- coding: utf-8 -*-
# preproc.py
# made by Min-Seok Kwon
# 2019-11-04 07:55:02
#########################
import sys
import os
import argparse
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path = "/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path = "/ms1/bin/python_lib"
else:
    sys_path = "/home/mk446/bin/python_lib"
sys.path.append(sys_path)

VERSION = "0.1"
VERSION_DATE = "2019.11.04"
PROG = "preproc"


def get_options():
    parser = argparse.ArgumentParser(
        usage='%(prog)s <sub-command> [options]', description='%(prog)s ver' + VERSION + " (" + VERSION_DATE + ")" + ': convert bam to image')
    parser.add_argument('-v', '--version', action='version',
                        version="%(prog)s ver" + VERSION + " (" + VERSION_DATE + ")")
    subparsers = parser.add_subparsers(title="sub-commands", dest="subcommand", metavar='', prog=PROG)

    p1 = subparsers.add_parser('preproc', help='quality check for VCF', description='quality check for VCF')
    p1.add_argument('-ds', dest='ds', default='', help='data source')
    p1.add_argument('-out', dest='out', default='', help='title of output file')
    p1.add_argument('-outtype', dest='outtype', default='', help='')
    # p1.add_argument('-ped', dest='ped', default='', help='title of output file')
    # p1.add_argument('-silence', dest='silence', action="store_true", default=False, help='don\'t print any log.')
    # p1.add_argument('-debug', dest='debug', action="store_true", default=False, help='turn on the debugging mode')

    # p1 = subparsers.add_parser('bamstat', help='quality metrics for BAM', description='quality metrics for BAM')
    # p1.add_argument('-bam', dest='bam', default='', help='BAM file')
    # p1.add_argument('-out', dest='out', default='', help='title of output file')
    # p1.add_argument('-silence', dest='silence', action="store_true", default=False, help='don\'t print any log.')
    # p1.add_argument('-debug', dest='debug', action="store_true", default=False, help='turn on the debugging mode')

    # p1 = subparsers.add_parser('vcfstat', help='quality metrics for VCF', description='quality metrics for VCF')
    # p1.add_argument('-vcf', dest='vcflist', default=[], help='VCF files', nargs="*")
    # p1.add_argument('-out', dest='out', default='', help='title of output file')
    # p1.add_argument('-ped', dest='ped', default='', help='title of output file')
    # p1.add_argument('-silence', dest='silence', action="store_true", default=False, help='don\'t print any log.')
    # p1.add_argument('-debug', dest='debug', action="store_true", default=False, help='turn on the debugging mode')

    # p1 = subparsers.add_parser('famvcfstat', help='quality metrics for family VCF', description='quality metrics for VCF')
    # p1.add_argument('-vcf', dest='vcflist', default=[], help='VCF files', nargs="*")
    # p1.add_argument('-out', dest='out', default='', help='title of output file')
    # p1.add_argument('-ped', dest='ped', default='', help='title of output file')
    # p1.add_argument('-silence', dest='silence', action="store_true", default=False, help='don\'t print any log.')
    # p1.add_argument('-debug', dest='debug', action="store_true", default=False, help='turn on the debugging mode')

    if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1][0] != '-'):
        sys.argv.append('-h')
    opt = vars(parser.parse_args())
    return opt


def cli():
    opt = get_options()
    # print (opt)
    if opt['subcommand'] == 'annot' and 'vcf' in opt.keys() and opt['vcf'] != "":
        av = annot.AnnotVCF(opt)
        av.run()
    if opt['subcommand'] == 'vcfqc' and 'vcflist' in opt.keys() and opt['vcflist'] != []:
        qcb = vcfqc.QCBoardVCF(opt)
        qcb.run()
    if opt['subcommand'] == 'vcfstat' and 'vcflist' in opt.keys() and opt['vcflist'] != []:
        if opt['ped'] == "":
            for vcf in opt['vcflist']:
                opt['vcf'] = vcf
                qcb = vcfstat.QCBoardVCFSTAT(opt)
                qcb.run()
        else:
            qcb = famvcfstat.QCBoardFamilyVCFSTAT(opt)
            qcb.run()


if __name__ == "__main__":
    cli()


if __name__ == "__main__":
    cli()
