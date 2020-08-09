CGAP Example
============

Micro annotation
----------------

.. code::

    mutanno annot \
            -vcf input.vcf \
            -ds datastructure_microannot_v0.4.4.json \
            -out output.annot.vcf \
            -sourcefile microannot_datasource.v0.4.4_200614.tsi.gz \
            -split_multi_allelic_variant \
            -genoinfo \
            -single_source_mode

* ds file: s3://maestro-resources/MICROANNOT/datastructure_microannot_v0.4.4.json
* mutanno: https://github.com/dbmi-bgm/mutanno/releases/tag/0.4.1 


Full annotation
---------------

.. code::

    mutanno annot \
            -vcf input.vcf \
            -ds datastructure_fullannot_v0.4.6.json \
            -out output.vcf \
            -sourcefile fullannot_source_file.mti.gz \
            -hg19 \
            -chain hg38ToHg19.over.chain.gz \
            -clean_tag MUTANNO SpliceAI CLINVAR gnomADgenome

* ds file: 
* mutanno: https://github.com/dbmi-bgm/mutanno/releases/tag/0.4.2
* source file: s3://maestro-resources/FULLANNO/merged.mti.gz

Gene annotation
---------------

.. code::

    mutanno makedata \
            -ds tests/data/datastructure_gene_v0.4.6ds.json \
            -out mvp_gene_datasource_v0.4.6.coding_gene_main_chrom \
            -vartype CODING_GENE_MAIN_CHROM \
            -outtype json

    gzip -c mvp_gene_datasource_v0.4.6.coding_gene_main_chrom.json > mvp_gene_datasource_v0.4.6.coding_gene_main_chrom.json.gz

* data file: https://www.dropbox.com/s/s6ahfq0gdn99uu8/mvp_gene_datasource_v0.4.6.coding_gene_main_chrom.json.gz?dl=0