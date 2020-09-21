Run Command in CGAP Project
===========================

Convert VEP to .mti for novel InDels
------------------------------------

Run VEP
^^^^^^^

.. code::

    ./bin/ensembl-vep-release-99/vep \
            -i input_novel_indels.vcf \
            -o input_novel_indels.vep.vcf \ 
            --hgvs \
            --fasta GRCh38_full_analysis_set_plus_decoy_hla.fa \
            --assembly GRCh38 \
            --use_given_ref \
            --offline \
            --cache_version 99 \
            --dir_cache ./bin/nonindexed_vep_cache/homo_sapiens_vep \
            --everything \
            --force_overwrite \
            --vcf \
            --plugin MaxEntScan,./bin/VEP_plugins-release-99/fordownload \
            --plugin TSSDistance \
            --dir_plugins ./bin/VEP_plugins-release-99 \
            --plugin SpliceRegion,Extended


.. note:: 
   
    `How to install VEP. <https://uswest.ensembl.org/info/docs/tools/vep/script/vep_download.html>`_

    1. download `VEP file (v99) from ENSEMBL <ftp://ftp.ensembl.org/pub/release-99/variation/vep/homo_sapiens_vep_99_GRCh38.tar.gz>`_ .
    2. untar and ungzip the downloaded file. 
    3. install plugins (https://uswest.ensembl.org/info/docs/tools/vep/script/vep_plugins.html)


For the micro annotation
^^^^^^^^^^^^^^^^^^^^^^^^

.. code::

    mutanno preprocess \
            -infile input_novel_indels.vep.vcf \
            -ds datastructure_microannot_v0.4.4.json \
            -out additional_novel_indels.vep.microannot.mti \
            -vep2mti

    bgzip -c additional_novel_indels.vep.microannot.mti > additional_novel_indels.vep.microannot.mti.gz
    tabix -f -p vcf additional_novel_indels.vep.microannot.mti.gz;

For the full annotation
^^^^^^^^^^^^^^^^^^^^^^^

.. code::

    mutanno preprocess \
            -infile input_novel_indels.vep.vcf \
            -out additional_novel_indels.vep.fullannot.mti \
            -vep2mti
    
    bgzip -c additional_novel_indels.vep.fullannot.mti > additional_novel_indels.vep.fullannot.mti.gz
    tabix -f -p vcf additional_novel_indels.vep.fullannot.mti.gz;

Micro annotation
----------------

.. code::

    mutanno annot \
            -vcf input.vcf \
            -ds datastructure_microannot_v0.4.5.json \
            -out output.annot.vcf \
            -sourcefile microannot_datasource.v0.4.4_200614.mti.gz additional.mti.gz \
                additional_novel_indels.vep.microannot.mti.gz \
            -split_multi_allelic_variant \
            -genoinfo \
            -single_source_mode

* ds file: https://github.com/dbmi-bgm/mutanno/blob/master/tests/data_structure_json/datastructure_microannot_v0.4.5.json
* mutanno: https://github.com/dbmi-bgm/mutanno/releases/tag/0.4.1 
* source file: 
    * s3://maestro-resources/MICROANNOT/microannot_datasource.v0.4.4_200614.tsi.gz and its tabix index file (.tbi)
    * Dropbox: 

.. tabs::

    .. tab:: input raw vcf
        
        .. code-block::
           :linenos:

            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12877_sample	NA12878_sample	NA12879_sample
            chr2	55544025	rs1045910	A	G	3047.94	.	AC=4;AF=0.667;AN=6;BaseQRankSum=0.502;DB;DP=148;ExcessHet=3.01;FS=1.374;MLEAC=4;MLEAF=0.667;MQ=60.00;MQRankSum=0.00;QD=20.59;ReadPosRankSum=0.549;SOR=0.709	GT:AD:DP:GQ:PL	0/1:27,20:47:99:534,0,756	1/1:0,50:50:99:1717,150,0	0/1:23,28:51:99:810,0,621

    .. tab:: output annotated vcf

        .. code-block::
           :linenos:

            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12877_sample	NA12878_sample	NA12879_sample
            chr2	55544025	rs1045910	A	G	3047.94	.	AC=4;AF=0.667;AN=6;BaseQRankSum=0.502;DB;DP=148;ExcessHet=3.01;FS=1.374;MLEAC=4;MLEAF=0.667;MQ=60.00;MQRankSum=0.00;QD=20.59;ReadPosRankSum=0.549;SOR=0.709;SAMPLEGENO=0/1|A/G|27/20|NA12877_sample,1/1|G/G|0/50|NA12878_sample,0/1|A/G|23/28|NA12879_sample;VEP=ENSG00000163001|ENST00000339012|Transcript|missense_variant|CFAP36|protein_coding,ENSG00000163001|ENST00000349456|Transcript|missense_variant|CFAP36|protein_coding,ENSG00000163001|ENST00000406691|Transcript|downstream_gene_variant|CFAP36|protein_coding,ENSG00000163001|ENST00000407816|Transcript|missense_variant~splice_region_variant|CFAP36|protein_coding,ENSG00000163001|ENST00000481791|Transcript|non_coding_transcript_exon_variant|CFAP36|retained_intron,ENSG00000163001|ENST00000490934|Transcript|non_coding_transcript_exon_variant|CFAP36|processed_transcript,ENSG00000275052|ENST00000611717|Transcript|downstream_gene_variant|PPP4R3B|protein_coding,ENSG00000275052|ENST00000616288|Transcript|downstream_gene_variant|PPP4R3B|protein_coding,ENSG00000275052|ENST00000616407|Transcript|downstream_gene_variant|PPP4R3B|protein_coding;gnomADgenome=9.40488e-01;SpliceAI=0.10	GT:AD:DP:GQ:PL	0/1:27,20:47:99:534,0,756	1/1:0,50:50:99:1717,150,0	0/1:23,28:51:99:810,0,621


Full annotation
---------------

.. code::

    mutanno annot \
            -vcf input.vcf \
            -ds datastructure_fullannot_v0.4.8.json \
            -out output.vcf \
            -sourcefile fullannot_source_file.mti.gz \
                additional_novel_indels.vep.fullannot.mti.gz \
            -hg19 \
            -chain hg38ToHg19.over.chain.gz \
            -clean_tag MUTANNO SpliceAI CLINVAR gnomADgenome

* ds file: https://github.com/dbmi-bgm/mutanno/blob/master/tests/data_structure_json/datastructure_fullannot_v0.4.8.json
* mutanno: https://github.com/dbmi-bgm/mutanno/releases/tag/0.4.3
* source file: s3://maestro-resources/FULLANNO/merged.mti.gz

Gene annotation
---------------

.. code::

    mutanno download \
            -source_path datasource_directory \
            -source all \
            -version lastest \
            -refversion hg38 \
            -websource mutanno

    mutanno makedata \
            -ds tests/data/datastructure_gene_v0.4.6ds.json \
            -out mvp_gene_datasource_v0.4.6.coding_gene_main_chrom \
            -vartype CODING_GENE_MAIN_CHROM \
            -outtype json

    gzip -c mvp_gene_datasource_v0.4.6.coding_gene_main_chrom.json > mvp_gene_datasource_v0.4.6.coding_gene_main_chrom.json.gz

* ds file: https://github.com/dbmi-bgm/mutanno/blob/master/tests/data_structure_json/datastructure_gene_v0.4.6ds.json
* out file: https://www.dropbox.com/s/s6ahfq0gdn99uu8/mvp_gene_datasource_v0.4.6.coding_gene_main_chrom.json.gz?dl=0