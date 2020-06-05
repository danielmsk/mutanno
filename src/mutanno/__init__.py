# -*- coding: utf-8 -*-
from ._logging import get_logger
from ._options import get_options
from .annot2 import annotvcf
from . import makegenedata
from . import makedata
from . import convert
from . import precal
from . import preprocess
from . import web
from . import validate


def cli():
    global _log

    opt = get_options()
    opt.log = get_logger(silence=opt.silence, debug=opt.debug, logfile=opt.logfile)
    dispatch_job(opt)


def dispatch_job(opt):
    if opt.subcommand == 'annot' and opt.vcf != "":
        an = annotvcf.VCFAnnotator(opt)
        an.run()
        
    if opt.subcommand == 'makedata' and  opt.out != "":
        opt.vartype = opt.vartype.upper()
        if 'GENE' in opt.vartype:
            md = makegenedata.GeneDataSourceFile(opt)
            md.make_single_source_file()
        else:
            md = makedata.DataSourceGenerator(opt)
            if opt.check:
                md.check_datasourcefile()
            else:
                if opt.blocksize == '':
                    md.make_single_source_file()
                else:
                    md.make_single_source_file_with_block()
    if opt.subcommand == 'convert' and opt.infile != "":
        if opt.ds != "":
            cv = convert.CONVFILE(opt)
            cv.convert()
        if opt.vcf2tsv:
            cv = convert.VCF2TSV(opt)
            cv.convert()
        if opt.vep2tab:
            cv = convert.VEP2TAB(opt)
            cv.convert()
    if opt.subcommand == 'validate':
        if opt.vcf != "":
            va = validate.VCFValidator(opt)
            va.validate()

    if opt.subcommand == 'precal':
        mp = precal.PreCalculate(opt)
        mp.run()
    if opt.subcommand == 'preprocess':
        if opt.make_dbnsfp_transcript:
            obj = preprocess.MakeDbnsfpTranscript(opt.infile, opt.out, opt.ds, opt)
            obj.run()
    if opt.subcommand == 'web' and opt.ds != "":
        web.runserver(opt)
