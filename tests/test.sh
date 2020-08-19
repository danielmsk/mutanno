# python ../src/mutanno.py annot \
# mutanno annot \
#         -vcf data/GAPFI3JX5D2J.vcf.gz \
#         -ds data_structure_json/archive/datastructure_fullannot_v0.4.6.json \
#         -out out/GAPFI3JX5D2J_fullannot_v0.4.6.vcf \
#         -sourcefile data/GAPFI3JX5D2J.vcf.gz_fullannot_source_v0.4.6.tsi.gz \
#         -hg19 \
#         -chain data/hg38ToHg19.over.chain.gz \
#         -clean_tag MUTANNO SpliceAI CLINVAR gnomADgenome
# vcf-validator out/GAPFI3JX5D2J_fullannot_v0.4.6.vcf

# python ../src/mutanno.py annot \
mutanno annot \
        -vcf data/GAPFI3JX5D2J.vcf.gz \
        -ds data_structure_json/datastructure_fullannot_v0.4.8.json \
        -out out/GAPFI3JX5D2J_fullannot_v0.4.8.vcf \
        -sourcefile data/GAPFI3JX5D2J.vcf.gz_fullannot_source_v0.4.8.tsi.gz \
        -hg19 \
        -chain data/hg38ToHg19.over.chain.gz \
        -clean_tag MUTANNO SpliceAI CLINVAR gnomADgenome
vcf-validator out/GAPFI3JX5D2J_fullannot_v0.4.8.vcf
python validate_vcf.py out/GAPFI3JX5D2J_fullannot_v0.4.8.vcf

## CLINVAR example
# python ../src/mutanno.py annot \
#         -vcf data/GAPFI3JX5D2J_CLINVAR.vcf.gz \
#         -ds data_structure_json/datastructure_fullannot_v0.4.8.json \
#         -out out/GAPFI3JX5D2J_CLINVAR_fullannot_v0.4.8.vcf \
#         -sourcefile data/GAPFI3JX5D2J.vcf.gz_fullannot_source_v0.4.8.tsi.gz \
#         -hg19 \
#         -chain data/hg38ToHg19.over.chain.gz \
#         -clean_tag MUTANNO SpliceAI CLINVAR gnomADgenome

## CLINVAR only 1 example
# python ../src/mutanno.py annot \
#         -vcf data/GAPFI3JX5D2J_CLINVAR_1.vcf.gz \
#         -ds data_structure_json/datastructure_fullannot_v0.4.8.json \
#         -out out/GAPFI3JX5D2J_CLINVAR_1_fullannot_v0.4.8.vcf \
#         -sourcefile data/GAPFI3JX5D2J.vcf.gz_fullannot_source_v0.4.8.tsi.gz \
#         -hg19 \
#         -chain data/hg38ToHg19.over.chain.gz \
#         -clean_tag MUTANNO SpliceAI CLINVAR gnomADgenome


## Multi-allelic variant
# python ../src/mutanno.py annot \
#         -vcf data/GAPFI3JX5D2J_multiallele.vcf.gz \
#         -ds data_structure_json/datastructure_fullannot_v0.4.8.json \
#         -out out/GAPFI3JX5D2J_multiallele_fullannot_v0.4.8.vcf \
#         -sourcefile data/GAPFI3JX5D2J.vcf.gz_fullannot_source_v0.4.8.tsi.gz \
#         -hg19 \
#         -chain data/hg38ToHg19.over.chain.gz \
#         -clean_tag MUTANNO SpliceAI CLINVAR gnomADgenome
# vcf-validator out/GAPFI3JX5D2J_multiallele_fullannot_v0.4.8.vcf

## Multi-allelic variant for 3 variants
# python ../src/mutanno.py annot \
#         -vcf data/GAPFI3JX5D2J_multiallele_1.vcf.gz \
#         -ds data_structure_json/datastructure_fullannot_v0.4.8.json \
#         -out out/GAPFI3JX5D2J_multiallele_1_fullannot_v0.4.8.vcf \
#         -sourcefile data/GAPFI3JX5D2J.vcf.gz_fullannot_source_v0.4.8.tsi.gz \
#         -hg19 \
#         -chain data/hg38ToHg19.over.chain.gz \
#         -clean_tag MUTANNO SpliceAI CLINVAR gnomADgenome
# vcf-validator out/GAPFI3JX5D2J_multiallele_1_fullannot_v0.4.8.vcf

# python ../src/mutanno.py annot \
#         -vcf data/test_trio_multiallele.vcf.gz \
#         -ds data_structure_json/datastructure_fullannot_v0.4.8.json \
#         -out out/test_trio_multiallele_fullannot_v0.4.8.vcf \
#         -sourcefile data/GAPFI3JX5D2J.vcf.gz_fullannot_source_v0.4.8.tsi.gz \
#         -hg19 \
#         -chain data/hg38ToHg19.over.chain.gz \
#         -clean_tag MUTANNO SpliceAI CLINVAR gnomADgenome



## Micro-annotation

# mutanno annot \
#         -vcf data/test_trio_1000_multiallele.vcf.gz \
#         -ds data_structure_json/datastructure_microannot_v0.4.4.json \
#         -out out/test_trio_multiallele_microannot_v0.4.8.vcf \
#         -sourcefile data/test_trio_multiallele.vcf.gz_microannot_datasource.v0.4.4_200614.tsi.gz \
#         -split_multi_allelic_variant \
#         -genoinfo \
#         -single_source_mode

# mutanno annot \
#         -vcf data/test_trio_1000_multiallele_2.vcf.gz \
#         -ds data_structure_json/datastructure_microannot_v0.4.4.json \
#         -out out/test_trio_multiallele_microannot_v0.4.4.vcf \
#         -sourcefile data/test_trio_1000_multiallele_2.vcf.gz_microannot_source_v0.4.4.tsi.gz \
#         -split_multi_allelic_variant \
#         -genoinfo \
#         -single_source_mode
