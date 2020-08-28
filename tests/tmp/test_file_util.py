#!/usr/bin/env python
# -*- coding: utf-8 -*-
from mutanno.util import file_util


def test_zero_format():
    assert file_util.zero_format(12, 5) == "00012"
    assert file_util.zero_format(1234, 2) == "1234"


# def test_load_json():
#     testjson = "datastructure.json"
#     assert file_util.load_json(testjson) == ""


def test_decodeb():
    assert file_util.decodeb("abcd") == "abcd"
    assert file_util.decodeb(b"abcd") == "abcd"


def test_rmext():
    assert file_util.rmext("abcd.vcf.gz", ".vcf.gz") == "abcd"
