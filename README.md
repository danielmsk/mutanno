# MutAnno
[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![Build Status](https://travis-ci.org/dbmi-bgm/mutanno.svg?branch=master)](https://travis-ci.org/dbmi-bgm/mutanno)


* github : https://github.com/dbmi-bgm/mutanno
* manual : https://mutanno.readthedocs.io/en/latest/

```
mutanno annot -vcf trio_test2.vcf \
    -out trio_test2.annot.vcf \
    -ds tests/datastructure_microannot_v1.0.json \
    -sourcefile sourcefile.tsv.gz \
    -split_multi_allelic_variant \
    -blocksize 10000
```

