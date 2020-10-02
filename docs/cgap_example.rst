Run Command in CGAP Project
===========================

.. toctree::
   :numbered:
   :maxdepth: 4

   cgap_example


Prerequisites
-------------

* vep (v99)
* tabix
* bgzip
* vcf-sort (vcf-tools)



Make Single Data Source File
----------------------------

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


Download and Preprocess Source for Micro Annotation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

VEP
***

* run VEP first.

.. code::

    mutanno preprocess \
            -infile input.vep.vcf \
            -ds datastructure_microannot_v0.4.4.json \
            -out vep.microannot.mti \
            -vep2mti

    bgzip -c additional_novel_indels.vep.microannot.mti > additional_novel_indels.vep.microannot.mti.gz
    tabix -f -p vcf additional_novel_indels.vep.microannot.mti.gz;


gnomAD
******

* input files
   * `gnomad.genomes.r3.0.sites.chr#CHROM#.vcf.bgz` and their .tbi file
* output file
   * `CLINVAR_hg38_20200927_variant.mti.gz` and its .tbi file
   * `CLINVAR_hg38_20200927_submission.sorted.mti.gz` and its .tbi file
* execution time: ~2hr

.. code::

    mutanno download \
            -source_path datasource_directory \
            -source gnomad \
            -version latest \
            -refversion hg38




ClinVar
*******

* input files
   * `clinvar_20200927.vcf.gz` (32 MB)
   * `ClinVarFullRelease_2020-09.xml.gz` (1.1 GB)
* output file
   * `CLINVAR_hg38_20200927_variant.mti.gz` and its .tbi file
   * `CLINVAR_hg38_20200927_submission.sorted.mti.gz` and its .tbi file
* execution time: ~1hr

.. code::

    mutanno download \
            -source_path datasource_directory \
            -source clinvar \
            -version latest \
            -refversion hg38

SpliceAI
********

* Download SpliceAI file from Illumina site
* input files
   * `spliceai_scores.raw.snv.hg38.vcf.gz` (28.8 GB)
   * `spliceai_scores.raw.snv.hg38.vcf.gz.tbi`
* output file
   * `SPLICEAI_hg38.mti.gz`
* execution time: ~1hr

.. code::

    mutanno preprocess \
            -infile datasource_directory/download/SPLICEAI/hg38/spliceai_scores.raw.snv.hg38.vcf.gz \
            -ds mutanno/tests/data_structure_json/datastructure_microannot_v0.4.5.json \
            -out datasource_directory/SPLICEAI_hg38.mti \
            -spliceai2mti


Make single source
******************

.. code::

    mutanno makedata \
            -ds mutanno/tests/data_structure_json/datastructure_microannot_v0.4.5ds.json \
            -out microannot_datasource.tsi \
            -vartype SNV \
            -blocksize 10000 \
            -region 1:1-100000 


Download and Preprocess Source for Full Annotation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


CYTOBAND
********

* input files
   * `spliceai_scores.raw.snv.hg38.vcf.gz` (28.8 GB)
   * `spliceai_scores.raw.snv.hg38.vcf.gz.tbi`
* output file
   * `SPLICEAI_hg38.mti.gz`
* execution time: ~1min

.. code::

    mutanno download \
            -source_path datasource_directory \
            -source CYTOBAND \
            -version latest \
            -refversion hg38


GENOMIC_SUPER_DUPLICATES
************************

* input files
   * `genomicSuperDups.txt.gz` (4 MB)
   * `genomicSuperDups.sql` (2 KB)
* output file
   * `GENOMIC_SUPER_DUPLICATES_hg38_2014-10-19.bed.gz` (4 MB)
* execution time: ~1min

.. code::

    mutanno download \
            -source_path datasource_directory \
            -source GENOMIC_SUPER_DUPLICATES \
            -version latest \
            -refversion hg38


SIMPLE_REPEAT
*************

* input files
   * `genomicSuperDups.txt.gz` (4 MB)
   * `genomicSuperDups.sql` (2 KB)
* output file
   * `SIMPLE_REPEAT_hg38_2019-03-11.bed.gz` (4 MB)
* execution time: ~1min

.. code::

    mutanno download \
            -source_path datasource_directory \
            -source SIMPLE_REPEAT \
            -version latest \
            -refversion hg38

RMSK
****

* input files
   * `genomicSuperDups.txt.gz` (4 MB)
   * `genomicSuperDups.sql` (2 KB)
* output file
   * `RMSK_hg38_2019-03-11.bed.gz` (165 MB)
* execution time: ~3min

.. code::

    mutanno download \
            -source_path datasource_directory \
            -source RMSK \
            -version latest \
            -refversion hg38


NESTED_REPEATS
**************

* input files
   * `genomicSuperDups.txt.gz` (4 MB)
   * `genomicSuperDups.sql` (2 KB)
* output file
   * `GENOMIC_SUPER_DUPLICATES_hg38_2014-10-19.bed.gz` (4 MB)
* execution time: ~1min

.. code::

    mutanno download \
            -source_path datasource_directory \
            -source NESTED_REPEATS \
            -version latest \
            -refversion hg38


MICROSATELLITE
**************

* input files
   * `genomicSuperDups.txt.gz` (4 MB)
   * `genomicSuperDups.sql` (2 KB)
* output file
   * `GENOMIC_SUPER_DUPLICATES_hg38_2014-10-19.bed.gz` (4 MB)
* execution time: ~1min

.. code::

    mutanno download \
            -source_path datasource_directory \
            -source MICROSATELLITE \
            -version latest \
            -refversion hg38


UNIPROT_TRANSMEM
****************

* input files
   * `genomicSuperDups.txt.gz` (4 MB)
   * `genomicSuperDups.sql` (2 KB)
* output file
   * `UNIPROT_TRANSMEM_hg38_2020-08-20.bed.gz` (4 MB)
* execution time: ~1min

.. code::

    mutanno download \
            -source_path datasource_directory \
            -source UNIPROT_TRANSMEM \
            -version latest \
            -refversion hg38


DBSNP
*****

* input files
   * `GRCh38_latest_dbSNP_all.vcf.gz` (16 GB)
* output file
   * `UNIPROT_TRANSMEM_hg38_2020-08-20.bed.gz` (4 MB)
* execution time: ~1hr

.. code::

    mutanno download \
            -source_path datasource_directory \
            -source DBSNP \
            -version latest \
            -refversion hg38


PRIMATEAI
*********

* Download PrimateAI file from Illumina site
* input files
   * `PrimateAI_scores_v0.2_hg38.tsv.gz` (868 MB)
* output file
   * `SPLICEAI_hg38.mti.gz`
* execution time: ~1hr

.. code::

    mutanno preprocess \
            -infile datasource_directory/download/PRIMATEAI/hg38/PrimateAI_scores_v0.2_hg38.tsv.gz \
            -ds mutanno/tests/data_structure_json/datastructure_microannot_v0.4.5.json \
            -out datasource_directory/PrimateAI_hg38.mti \
            -primateai2mti

TOPMED
******

* Download bravo-dbsnp-all.tsv.gz (hg19) file from Bravo site (https://bravo.sph.umich.edu)
   * There is no hg38 version for topmed. 
   * In this preprocessing step, we can lift over to hg38 from hg19.
* input files
   * `bravo-dbsnp-all.tsv.gz` (6.1 GB)
* output file
   * `TOPMED_hg38.mti.gz`
* execution time: ~1hr
   * liftover
   * vcf-sort
   * tabixgz

.. code::

    mutanno preprocess \
            -infile datasource_directory/download/TOPMED/hg19/bravo-dbsnp-all.tsv.gz \
            -ds mutanno/tests/data_structure_json/datastructure_microannot_v0.4.5.json \
            -out datasource_directory/TOPMED_hg38.mti \
            -topmed2mti

CADD
****

* input files
   * `GRCh38_latest_dbSNP_all.vcf.gz` (16 GB)
* output file
   * `UNIPROT_TRANSMEM_hg38_2020-08-20.bed.gz` (4 MB)
* execution time: ~1hr

.. code::

    mutanno download \
            -source_path datasource_directory \
            -source DBSNP \
            -version latest \
            -refversion hg38


UK10K
*****
* input files
   * `GRCh38_latest_dbSNP_all.vcf.gz` (16 GB)
* output file
   * `UNIPROT_TRANSMEM_hg38_2020-08-20.bed.gz` (4 MB)
* execution time: ~1hr

.. code::

    mutanno download \
            -source_path datasource_directory \
            -source UK10K \
            -version latest \
            -refversion hg38


COSMIC
******

* Download CosmicMutantExport.tsv.gz file from Cosmic site (https://cancer.sanger.ac.uk/cosmic/download; login required)
* input files
   * `CosmicMutantExport.tsv.gz` (6.5 GB)
* output file
   * `SPLICEAI_hg38.mti.gz`
* execution time: ~1hr

.. code::

    mutanno preprocess \
            -infile datasource_directory/download/TOPMED/hg38/bravo-dbsnp-all.vcf.gz \
            -ds mutanno/tests/data_structure_json/datastructure_microannot_v0.4.5.json \
            -out datasource_directory/TOPMED_hg38.mti \
            -cosmic2mti



CONSERVATION
************

* input files
   * `GRCh38_latest_dbSNP_all.vcf.gz` (16 GB)
* output file
   * `UNIPROT_TRANSMEM_hg38_2020-08-20.bed.gz` (4 MB)
* execution time: ~1hr

.. code::

    mutanno download \
            -source_path datasource_directory \
            -source PHASTCONS \
            -version 20way \
            -refversion hg38
    
    mutanno download \
            -source_path datasource_directory \
            -source PHASTCONS \
            -version 30way \
            -refversion hg38

    mutanno download \
            -source_path datasource_directory \
            -source PHASTCONS \
            -version 100way \
            -refversion hg38

    mutanno download \
            -source_path datasource_directory \
            -source PHYLOP \
            -version 20way \
            -refversion hg38
    
    mutanno download \
            -source_path datasource_directory \
            -source PHYLOP \
            -version 30way \
            -refversion hg38

    mutanno download \
            -source_path datasource_directory \
            -source PHYLOP \
            -version 100way \
            -refversion hg38

    mutanno download \
            -source_path datasource_directory \
            -source GERP \
            -version latest \
            -refversion hg38

    mutanno download \
            -source_path datasource_directory \
            -source SHIPHY \
            -version latest \
            -refversion hg38


    mutanno preprocess \
            -infile datasource_directory/PHASTCONS_hg38_20way.tsi.gz \
                    datasource_directory/PHASTCONS_hg38_30way.tsi.gz \
                    datasource_directory/PHASTCONS_hg38_100way.tsi.gz \
                    datasource_directory/PHYLOP_hg38_20way.tsi.gz \
                    datasource_directory/PHYLOP_hg38_30way.tsi.gz \
                    datasource_directory/PHYLOP_hg38_100way.tsi.gz \
                    datasource_directory/GERP_hg38.tsi.gz \
                    datasource_directory/SHIPHY_hg38.tsi.gz \
            -ds mutanno/tests/data_structure_json/datastructure_microannot_v0.4.5.json \
            -out datasource_directory/CONSERVATION_hg38.chr___CHROM___.bed.gz \
            -conservation2mti


MAX_POP_AF


Make single source
******************

.. code::

    mutanno makedata \
            -ds mutanno/tests/data_structure_json/datastructure_fullannot_v0.4.8.json \
            -out single_datasource_file.tsi \
            -vartype SNV \
            -blocksize 1000




Download Data Source for Gene Annotation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code::

    mutanno download \
            -source_path datasource_directory \
            -source gnomAD \
            -version latest \
            -refversion hg38 \

    mutanno download \
            -source_path datasource_directory \
            -source CLINVAR \
            -version latest \
            -refversion hg38 \

    mutanno download \
            -source_path datasource_directory \
            -source UK10K \
            -version latest \
            -refversion hg38 \



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


Make Additional .mti for Novel InDels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code::

    mutanno preprocess \
            -infile input_novel_indels.vep.vcf \
            -ds datastructure_microannot_v0.4.4.json \
            -out additional_novel_indels.vep.microannot.mti \
            -vep2mti

    bgzip -c additional_novel_indels.vep.microannot.mti > additional_novel_indels.vep.microannot.mti.gz
    tabix -f -p vcf additional_novel_indels.vep.microannot.mti.gz;


Make Additional .mti for Novel InDels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code::

    mutanno preprocess \
            -infile input_novel_indels.vep.vcf \
            -out additional_novel_indels.vep.fullannot.mti \
            -vep2mti
    
    bgzip -c additional_novel_indels.vep.fullannot.mti > additional_novel_indels.vep.fullannot.mti.gz
    tabix -f -p vcf additional_novel_indels.vep.fullannot.mti.gz;



Annotation
----------


Run Micro-Annotation
^^^^^^^^^^^^^^^^^^^^


.. code::

    mutanno annot \
            -vcf input.vcf \
            -ds datastructure_microannot_v0.4.5.json \
            -out output.annot.vcf \
            -sourcefile microannot_datasource.v0.4.4_200614.mti.gz additional.mti.gz \
                additional_novel_indels.vep.microannot.mti.gz \
            -split_multi_allelic_variant \
            -genoinfo \
            -use_raw_source

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

Run Full-Annotation
^^^^^^^^^^^^^^^^^^^

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
            -version latest \
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



