# -*- coding: utf-8 -*-
from ._logging import get_logger
from ._options import get_options
from . import annotvcf
from .download import Downloader
from . import makegenedata
from . import makedata
from . import convert
from . import precal
from . import preprocess
from . import web
from . import validate
from . import viewer
from .model.datastructure import DataSourceListStructure


def cli():
    global _log

    opt = get_options()
    opt.log = get_logger(silence=opt.silence, debug=opt.debug, logfile=opt.logfile)
    dispatch_job(opt)


def update_option(opt):
    dslist = DataSourceListStructure(opt.ds)
    ds_field_list = list(dslist.__dict__.keys())
    for k1 in opt.__dict__.keys():
        if k1 in ds_field_list:
            if dslist.__dict__[k1] != "":
                opt.__dict__[k1] = dslist.__dict__[k1]
    return opt


def dispatch_job(opt):
    if opt.subcommand == 'annot' and opt.vcf != "":
        opt = update_option(opt)
        an = annotvcf.VCFAnnotator(opt)
        an.run()

    if opt.subcommand == 'download' and opt.dir != "":
        dn = Downloader(opt)
        dn.run()

    if opt.subcommand == 'makedata' and opt.out != "":
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

    if opt.subcommand == 'view':
        if opt.vcf != "":
            v1 = viewer.ANNOTVCFViewer(opt)
            v1.run()
