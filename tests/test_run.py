#!/usr/bin/env python
# -*- coding: utf-8 -*-
import test_conf
import mutanno
from mutanno.util import file_util


def test_run():
    opt = {
      'subcommand': 'annot', 
      'vcf': test_conf.TEST_VCF1,
      'ds': test_conf.TEST_DS1, 
      'sourcefile': test_conf.TEST_DSFILE1,
      'blocksize': 100, 
      'add_genoinfo': True, 
      'split_multi_allelic_variant': True,
      'load_source_in_memory': False, 
      'out': test_conf.TEST_OUT1, 
      'clean_tag_list': '',
      'remove_unannotated_variant': False
    }
    mutanno.dispatch_job(opt)
    assert file_util.fileOpen(opt['out']) == file_util.fileOpen(opt['out'] + ".0")
    # assert test_conf.TEST_VCF1 == ""
    # assert os.path.abspath('./') == "/"
    # assert test_conf.TESTS_DIR == ""
    pass
