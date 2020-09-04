Datasource List
===============

All preprocessed data is available in `dropbox shared folder <https://www.dropbox.com/sh/trjwttkf0ypn8c8/AAD3yVK-HBbkm0eXsYSq5r85a?dl=0>`_.


**Variant annotation (hg38)**

=========== =============  =========== ========== ========= ===============
Source name Category       Version     Date       Download  link
----------- -------------  ----------- ---------- --------- ---------------
hg38 (GRCh38)
---------------------------------------------------------------------------
VEP         Annotation     v99         11/04/2019 manually  `⬇ <#vep>`_
dbSNP       PopulationDB               00/00/0000 auto      `⬇ <#dbsnp>`_
gnomAD      PopulationDB   3.0                    auto      `⬇ <#gnomad>`_
UK10K       PopulationDB   20160215    02/15/2016 auto      `⬇ <#uk10k>`_
CLINVAR     VariantDB      20200106    01/06/2020 auto      `⬇ <#clinvar>`_
SpliceAI    Pathogenicity  20191004    01/06/2019 auto      `⬇ <#spliceai>`_
CADD        Pathogenicity  00000       00/00/0000 auto      `⬇ <#cadd>`_
ENSEMBL     Pathogenicity  00000       00/00/0000 auto      `⬇ <#ensembl>`_
GERP        Conservation   100_mammals 01/01/2000 auto      `⬇ <#phylop>`_
PHASTCONS   Conservation   100way      01/01/2000 auto      `⬇ <#phylop>`_
PHYLOP      Conservation   20way       04/16/2015 auto      `⬇ <#phylop>`_
PHYLOP      Conservation   30way       11/05/2017 auto      `⬇ <#phylop>`_
PHYLOP      Conservation   100way      05/07/2015 auto      `⬇ <#phylop>`_
SIPHY       Conservation   20way       01/01/2000 auto      `⬇ <#shiphy>`_

hg19 (GRCh37)
---------------------------------------------------------------------------


=========== =============  =========== ========== ========= ===============

**auto**: support to download and preprocess automatically in mutanno


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
    
      wget http://
      mutanno preprocess \
              -source_path datasource_directory \
              -source phylop \
              -version 20way \
              -refversion hg38 \
              -websource mutanno


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

   .. code-block::
      :linenos:
      
      vep OOOO

3. preprocess VEP result (convert .mti)

   .. code-block::
      :linenos:
      
      mutanno preprocess \
              -out test.vcf



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

1. download data file (.bw) from `ensembl ftp <ftp://ftp.ensembl.org/pub/current_compara/conservation_scores/100_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw>`_
2. Convert .bw file to .wig using `bigWigToWig <http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig>`_
3. run following preprocessing comamnd. (.wig -> .mti.gz)

   .. code:: console
    
      mutanno preprocess


.. note::

   The current GERP version is 111_mammals (veriosn date is 7/18/20/). This part needs to be updated.


PHASTCONS
^^^^^^^^^

1. download data file from `USCS database <ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr1.phastCons100way.wigFix.gz>`_
3. run following preprocessing comamnd. (.wig -> .mti.gz)

   .. code::
    
      mutanno preprocess

PHYLOP
^^^^^^

* web resource: `UCSC database phyloP100way <ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/>`_, `phyloP30way <ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP30way/>`_, `phyloP20way <ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP20way/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.

* download and preprocess automatically.

   .. code-block::
      :linenos:
      :emphasize-lines: 3,8,13
    
      mutanno download \
              -source_path datasource_directory \
              -source phylop \
              -version 20way \
              -refversion hg38
      mutanno download \
              -source_path datasource_directory \
              -source phylop \
              -version 30way \
              -refversion hg38
       mutanno download \
              -source_path datasource_directory \
              -source phylop \
              -version 100way \
              -refversion hg38   
              

* download preprocessed file from mutanno dropbox.

   .. code-block::
      :linenos:
      :emphasize-lines: 4,6,10,12,16,18
    
      mutanno download \
              -source_path datasource_directory \
              -source phylop \
              -version 20way \
              -refversion hg38 \
              -websource mutanno
      mutanno download \
              -source_path datasource_directory \
              -source phylop \
              -version 30way \
              -refversion hg38 \
              -websource mutanno
       mutanno download \
              -source_path datasource_directory \
              -source phylop \
              -version 100way \
              -refversion hg38 \
              -websource mutanno


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

* web resource: `gnomAD broser <https://gnomad.broadinstitute.org/downloads>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


SpliceAI
^^^^^^^^

* web resource: `gnomAD broser <https://gnomad.broadinstitute.org/downloads>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.



Variant database
----------------


CLINVAR
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
