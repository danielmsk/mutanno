#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import json
from mutanno.util import file_util
from mutanno.util import proc_util
import test_conf

def test_reverse_full_annotation():
    ds = test_conf.get_ds()
    dsjson = ds['fullannot']['ds']
    dslist = file_util.load_json(dsjson)
    for s1 in dslist['source']:
        print(s1['name'])




if __name__ == "__main__":
    test_reverse()

