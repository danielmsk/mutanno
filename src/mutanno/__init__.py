# -*- coding: utf-8 -*-
from ._logging import get_logger
from ._options import get_options
from . import annot
from . import makegenedata
from . import makedata
from . import convert
from . import precal
from . import preprocess
from . import web


def cli():
    global _log

    opt = get_options()
    opt['log'] = get_logger(silence=opt['silence'], debug=opt['debug'], logfile=opt['logfile'])
    dispatch_job(opt)


def dispatch_job(opt):
    if opt['subcommand'] == 'annot' and 'vcf' in opt.keys() and opt['vcf'] != "":
        annotvcf = annot.AnnotVCF(opt)
        if opt['load_source_in_memory']:
            annotvcf.run_by_loading_source()
        else:
            annotvcf.run()
    if opt['subcommand'] == 'makedata' and 'out' in opt.keys() and opt['out'] != "":
        opt['vartype'] = opt['vartype'].upper()
        if 'GENE' in opt['vartype']:
            md = makegenedata.GeneDataSourceFile(opt)
            md.make_single_source_file()
        else:
            md = makedata.DataSourceFile(opt)
            if opt['check']:
                md.check_datasourcefile()
            else:
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
    if opt['subcommand'] == 'precal':
        mp = precal.PreCalculate(opt)
        mp.run()
    if opt['subcommand'] == 'preprocess':
        if opt['make_dbnsfp_transcript']:
            obj = preprocess.MakeDbnsfpTranscript(opt['infile'], opt['out'], opt['ds'], opt)
            obj.run()
    if opt['subcommand'] == 'web' and 'ds' in opt.keys() and opt['ds'] != "":
        web.runserver(opt)
