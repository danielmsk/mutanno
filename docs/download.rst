Download Data Source (``download``)
===================================

To annotate variants, **MutAnno** need data source files and can download them with ``download`` sub-command.

.. code::

   mutanno download \
           -source_path datasource_directory \
           -source all \
           -version lastest \
           -refversion hg38 \
           -websource mutanno

* ``-source_path`` : data source path (default: `mutanno_source`)
* ``-source`` : source name (default: `all`)
* ``-version`` : source version (default: `lastest`)
* ``-refversion`` : reference version [`hg19`, `hg38`] (default: `hg38`)
* ``-websource`` : web source indication [``, `mutanno`] (default: ``). When `-websource mutanno` is used, only preprocessed files can be download from MutAnno web source.

.. note::

   If ``-source`` is ``all``, all data sources can be downloaded and preprocessed.


Datasource List
---------------

All preprocessed data is available in `dropbox shared folder <https://www.dropbox.com/sh/trjwttkf0ypn8c8/AAD3yVK-HBbkm0eXsYSq5r85a?dl=0>`_.


**Variant annotation (hg38)**

============================ =============  =========== ========== ============ ===============
**hg38 (GRCh38)**
-----------------------------------------------------------------------------------------------
**Source name**              **Category**   **Version** **Date**   **Download** **link**
---------------------------- -------------  ----------- ---------- ------------ ---------------
**VEP**                      Annotation     v101        08/26/2020 manually     `⬇ <#vep>`_
**VEP** *                    Annotation     v99         11/04/2019 manually     `⬇ <#vep>`_
**dbSNP** *                  PopulationDB               00/00/0000 auto         `⬇ <#dbsnp>`_
**gnomAD** *                 PopulationDB   3.0                    auto         `⬇ <#gnomad>`_
**UK10K** *                  PopulationDB   20160215    02/15/2016 auto         `⬇ <#uk10k>`_
**TOPMED** *                 PopulationDB   freeze 5    08/28/2017 auto         `⬇ <#topmed>`_
**CLINVAR** *                VariantDB      20200106    01/06/2020 auto         `⬇ <#clinvar>`_
**COSMIC** *                 VariantDB      v90         08/06/2019 auto         `⬇ <#cosmic>`_
**SPLICEAI** *               Pathogenicity  20191004    01/06/2019 manually     `⬇ <#spliceai>`_
**PRIMATEAI** *              Pathogenicity  v0.2_hg38   12/18/2019 manually     `⬇ <#primateai>`_
**CADD**                     Pathogenicity  1.6         03/26/2020 auto         `⬇ <#cadd>`_
**CADD** *                   Pathogenicity  1.5         02/22/2019 auto         `⬇ <#cadd>`_
**GERP** *                   Conservation   100_mammals 01/01/2000 auto         `⬇ <#phylop>`_
**PHASTCONS** *              Conservation   100way      07/17/2017 auto         `⬇ <#phylop>`_
**PHASTCONS** *              Conservation   30way       07/17/2017 auto         `⬇ <#phylop>`_
**PHASTCONS** *              Conservation   20way       07/17/2017 auto         `⬇ <#phylop>`_
**PHYLOP** *                 Conservation   100way      04/16/2015 auto         `⬇ <#phylop>`_
**PHYLOP** *                 Conservation   30way       11/05/2017 auto         `⬇ <#phylop>`_
**PHYLOP** *                 Conservation   20way       05/07/2015 auto         `⬇ <#phylop>`_
**SIPHY** *                  Conservation   20way       01/01/2000 auto         `⬇ <#shiphy>`_
**SUPER_DUPLICATES** *       Repeatitive    20way       01/01/2000 auto         `⬇ <#super_duplicates>`_
**SIMPLE_REPEAT** *          Repeatitive    20way       01/01/2000 auto         `⬇ <#simple_repeat>`_
**RMSK** *                   Repeatitive    20way       01/01/2000 auto         `⬇ <#rmsk>`_
**NESTED_REPEATS** *         Repeatitive    20way       01/01/2000 auto         `⬇ <#nested_repeats>`_
**MICROSATELLITE** *         Repeatitive    20way       08/23/2015 auto         `⬇ <#microsatellite>`_
============================ =============  =========== ========== ============ ===============

* **auto**: support to download and preprocess automatically in mutanno
* **star(*)**: the star(*) means 

=============== =============  =========== ========== ============ ===============
**hg19 (GRCh37)**
----------------------------------------------------------------------------------
**Source name** **Category**   **Version** **Date**   **Download** **link**
--------------- -------------  ----------- ---------- ------------ ---------------
**CADD**        Pathogenicity  1.6         03/26/2020 auto         `⬇ <#cadd>`_
=============== =============  =========== ========== ============ ===============


**Gene annotation**

=========== =============  =========== ========== ========= ===============
Source name Category       Version     Date       Download  link
----------- -------------  ----------- ---------- --------- ---------------
GTEx        Expression
GeneMetrics Conservation   20way       01/01/2000 auto      `⬇ <#shiphy>`_
=========== =============  =========== ========== ========= ===============

Download methods
----------------

1. download and preprocess automatically.

   .. code-block::
      :linenos:
      :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source phylop \
              -version 20way \
              -refversion hg38

2. download preprocessed file from mutanno dropbox

   `-websource mutanno` option doesn't run preprocessing module.

   .. code-block::
      :linenos:
      :emphasize-lines: 6
    
      mutanno download \
              -source_path datasource_directory \
              -source phylop \
              -version 20way \
              -refversion hg38 \
              -websource mutanno

3. download manually (using wget), and then run preprocess module.

   .. code-block::
      :linenos:
      :emphasize-lines: 6
    
      wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP20way/hg38.phyloP20way.bw

      mutanno preprocess \
              -infile datasource_directory/hg38.phyloP20way.bw \
              -ds phylop.datastructure.json \
              -out datasource_directory/hg38.phyloP20way.mti.gz


Variant annotation
------------------

VEP
^^^

* MutAnno doesn't support to download VEP raw data automatically, but support to download preprocessed files from MutAnno dropbox


Download preprocessed files from MutAnno dropbox
************************************************

   .. code-block::
      :linenos:
      :emphasize-lines: 3,6
    
      mutanno download \
              -source_path datasource_directory \
              -source vep \
              -version lastest \
              -refversion hg38 \
              -websource mutanno

Make VEP result files and then run preprocess
*********************************************

1. make mock vcf files

   .. code-block::
      :linenos:
      
      mutanno vcfmaker \
              -out test.vcf
      

2. run VEP


    * download ftp://ftp.ensembl.org/pub/release-99/variation/vep/homo_sapiens_vep_99_GRCh38.tar.gz

   .. code-block::
      :linenos:
      
      ./bin/ensembl-vep-release-99/vep \
            -i chr1_100001_200000.vcf \
            -o chr1_100001_200000.vcf.vep.txt \ 
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

3. preprocess VEP result (convert .mti)

   .. code-block::
      :linenos:
      
      mutanno preprocess \
              -infile datasource_directory/chr1_100001_200000.vcf.vep.txt \
              -ds vep.datastructure.json \
              -out datasource_directory/chr1_100001_200000.vcf.vep.mti.gz

   
   We can merge the chopped mti files into one single files using `vcf-merger`.



Population data
---------------

dbSNP
^^^^^

* web resource: `NCBI refseq <ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_dbSNP_all.vcf.gz>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


gnomAD
^^^^^^

* web resource: `gnomAD broser <https://gnomad.broadinstitute.org/downloads>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.
* For the hg19, v2.1.1 is available. And for the hg39, v3.0 is available.


UK10K
^^^^^
* web resource: `UK10K of Sanger institude <ftp://ngs.sanger.ac.uk/production/uk10k/UK10K_COHORT/REL-2012-06-02/UK10K_COHORT.20160215.sites.vcf.gz>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.
* Only hg19 version of UK10K is available. For the hg38 version, MutAnno do the liftover from hg19 in the preprocessing.


Conservation
------------

GERP
^^^^

* Download data file (.bw) from `ensembl ftp <ftp://ftp.ensembl.org/pub/current_compara/conservation_scores/100_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw>`_
* Convert .bw file to .wig using `bigWigToWig <http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


.. note::

   The current GERP version is 111_mammals (veriosn date is 7/18/20/). This part needs to be updated.


PHASTCONS
^^^^^^^^^

* Download data file from `USCS database <ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr1.phastCons100way.wigFix.gz>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.

PHYLOP
^^^^^^

* web resource: `UCSC database phyloP100way <ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/>`_, `phyloP30way <ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP30way/>`_, `phyloP20way <ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP20way/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


SIPHY
^^^^^

* web resource: `gnomAD broser <https://gnomad.broadinstitute.org/downloads>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


Pathogenicity
-------------

* web resource: `gnomAD broser <https://gnomad.broadinstitute.org/downloads>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


CADD
^^^^

* web resource: `CADD web source <https://krishna.gs.washington.edu/download/CADD/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


SpliceAI
^^^^^^^^

* web resource: `gnomAD broser <https://gnomad.broadinstitute.org/downloads>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.



Variant database
----------------


CLINVAR
^^^^^^^

* web resource: `NCBI ClinVar broser <https://www.ncbi.nlm.nih.gov/variation/docs/ClinVar_vcf_files/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.



COSMIC
^^^^^^^

* web resource: `gnomAD broser <https://gnomad.broadinstitute.org/downloads>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


https://www.ncbi.nlm.nih.gov/variation/docs/ClinVar_vcf_files/



Gene annotation
---------------


GTEx
^^^^

* web resource: `gnomAD broser <https://gnomad.broadinstitute.org/downloads>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.



ENSEMBL
^^^^^^^

* web resource: `gnomAD broser <https://gnomad.broadinstitute.org/downloads>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


