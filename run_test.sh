#!/bin/bash

PROG=mutanno
# PROG=python /home/mk446/bio/mutanno/SRC/src/mutanno.py

VCF=/home/mk446/bio/mutanno/TEST/GAPFI7HPIVFM.mod.vcf.gz
# OUT=/home/mk446/bio/mutanno/TEST/GAPFI7HPIVFM.mod.annot_v0.4.5.vcf
DS=/home/mk446/bio/mutanno/SRC/tests/data/datastructure_v0.4.5ds.json
SOURCEFILE=/home/mk446/bio/mutanno/DATASOURCE/MAINANNOT/mc_3k.tsi_mod.tsi.gz
CHAIN=/home/mk446/mutanno/TEST/hg38ToHg19.over.chain.gz

## MicroAnnotation
SOURCEFILE=/home/mk446/mutanno/DATASOURCE/MICROANNOT/v0.4.2_200512/microannot_datasource.v0.4.2_20200512.tsi.gz
OUT=/home/mk446/bio/mutanno/TEST/GAPFI7HPIVFM.mod.microannot.vcf
# $PROG \
python /home/mk446/bio/mutanno/SRC/src/mutanno.py \
    annot \
    -vcf $VCF \
    -out $OUT \
    -ds $DS \
    -sourcefile $SOURCEFILE \
    -genoinfo NA12877_sample NA12878_sample NA12879_sample \
    -hg19 \
    -chain $CHAIN \
    -split_multi_allelic_variant
