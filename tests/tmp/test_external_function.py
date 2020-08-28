#!/usr/bin/env python
# -*- coding: utf-8 -*-
# test_external_function.py
# made by Daniel Minseok Kwon
# 2020-06-01 13:22:37
#########################
import sys
sys.path.append('..')
from src.mutanno.util import file_util
from src.mutanno import external_functions
from src.mutanno.model.datasource import run_external_function



def test_external_function():
    function = "conv_consequence_delimiter"
    params = ["consequence"]
    paramvalues = {"consequence":"A&B"}
    assert run_external_function(function, params, paramvalues) == "A~B"
    paramvalues = {"consequence":"A-B"}
    assert run_external_function(function, params, paramvalues) == "A-B"

    function = "convert_vep_strand"
    params = ["vep_strand"]
    paramvalues = {"vep_strand":"-1"}
    assert run_external_function(function, params, paramvalues) == "0"

if __name__ == "__main__":
    test_external_function()
