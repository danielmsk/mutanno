# MutAnno
[![pypi version](https://img.shields.io/pypi/v/mutanno.svg)](https://pypi.org/project/mutanno/)
[![doc](https://readthedocs.org/projects/mutanno/badge/?version=latest)](https://mutanno.readthedocs.io/en/latest/)
<!-- [![pypi download](https://img.shields.io/pypi/dm/mutanno.svg)](https://pypi.org/project/mutanno/) -->
<!-- [![Build Status](https://travis-ci.org/dbmi-bgm/mutanno.svg?branch=master)](https://travis-ci.org/dbmi-bgm/mutanno) -->


* github : https://github.com/dbmi-bgm/mutanno
* manual : https://mutanno.readthedocs.io/en/latest/

For more details, see MutAnno [**documentation**](http://mutanno.readthedocs.io/en/latest).

```
mutanno annot -vcf trio_test2.vcf \
    -out trio_test2.annot.vcf \
    -ds tests/datastructure_microannot_v1.0.json \
    -sourcefile sourcefile.tsv.gz \
    -split_multi_allelic_variant \
    -blocksize 10000
```

