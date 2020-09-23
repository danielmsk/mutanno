Download Data Source (``download``)
===================================

To annotate variants, **MutAnno** need data source files and can download them with ``download`` sub-command.

.. code::

   mutanno download \
           -source_path datasource_directory \
           -source all \
           -version latest \
           -refversion hg38 \
           -websource mutanno

* ``-source_path`` : data source path (default: `mutanno_source`)
* ``-source`` : source name (default: `all`)
* ``-version`` : source version (default: `latest`)
* ``-refversion`` : reference version [`hg19`, `hg38`] (default: `hg38`)
* ``-websource`` : web source indication [``, `mutanno`] (default: ``). When `-websource mutanno` is used, only preprocessed files can be download from MutAnno web source.

.. note::

   If ``-source`` is ``all``, all data sources can be downloaded and preprocessed.


Datasource List
---------------

All preprocessed data is available in `dropbox shared folder <https://www.dropbox.com/sh/trjwttkf0ypn8c8/AAD3yVK-HBbkm0eXsYSq5r85a?dl=0>`_.


**Variant annotation (hg38)**

============================ ============= ============ ===============
**hg38 (GRCh38)**
-----------------------------------------------------------------------
**Source name**              **Category**  **Download** **link**
---------------------------- ------------- ------------ ---------------
**VEP** *                    Annotation    manually     `⬇ <#vep>`_
**dbSNP** *                  PopulationDB  auto         `⬇ <#dbsnp>`_
**gnomAD** *                 PopulationDB  auto         `⬇ <#gnomad>`_
**UK10K** *                  PopulationDB  auto         `⬇ <#uk10k>`_
**TOPMED** *                 PopulationDB  auto         `⬇ <#topmed>`_
**CLINVAR** *                VariantDB     auto         `⬇ <#clinvar>`_
**COSMIC** *                 VariantDB     auto         `⬇ <#cosmic>`_
**SPLICEAI** *               Pathogenicity manually     `⬇ <#spliceai>`_
**PRIMATEAI** *              Pathogenicity manually     `⬇ <#primateai>`_
**CADD** *                   Pathogenicity auto         `⬇ <#cadd>`_
**GERP** *                   Conservation  auto         `⬇ <#phylop>`_
**PHASTCONS** *              Conservation  auto         `⬇ <#phylop>`_
**PHYLOP** *                 Conservation  auto         `⬇ <#phylop>`_
**SIPHY** *                  Conservation  auto         `⬇ <#shiphy>`_
**SUPER_DUPLICATES** *       Repeatitive   auto         `⬇ <#super_duplicates>`_
**SIMPLE_REPEAT** *          Repeatitive   auto         `⬇ <#simple_repeat>`_
**RMSK** *                   Repeatitive   auto         `⬇ <#rmsk>`_
**NESTED_REPEATS** *         Repeatitive   auto         `⬇ <#nested_repeats>`_
**MICROSATELLITE** *         Repeatitive   auto         `⬇ <#microsatellite>`_
============================ ============= ============ ===============

* **auto**: support to download and preprocess automatically in mutanno
* **star(*)**: the star(*) means this version is used in CGAP project.

=============== ============= ============ ===============
**hg19 (GRCh37)**
----------------------------------------------------------
**Source name** **Category**  **Download** **link**
--------------- ------------- ------------ ---------------
**CADD**        Pathogenicity auto         `⬇ <#cadd>`_
=============== ============= ============ ===============


**Gene annotation**

============================== ============ ===============
**Source name**                **Download** **link**
------------------------------ ------------ ---------------
**ENSEMBLgene** *              auto         `⬇ <#ensemblgene>`_
**ENSEMBLgeneGRCh37** *        auto         `⬇ <#ensemblgenegrch37>`_
**CYTOBAND** *                 auto         `⬇ <#cytoband>`_
**RefSeq** *                   auto         `⬇ <#refseq>`_
**HGNC** *                     auto         `⬇ <#hgnc>`_
**ClinGen** *                  auto         `⬇ <#clingen>`_
**ClinGenDisease** *           auto         `⬇ <#clingendisease>`_
**ENSEMBLIDxrefTrscriptID** *  auto         `⬇ <#ensemblidxreftrscriptid>`_
**ENSEMBLIDxref** *            auto         `⬇ <#ensemblidxref>`_
**dbNSFP** *                   auto         `⬇ <#dbnsfp>`_
**gnomADmetrics** *            auto         `⬇ <#gnomadmetrics>`_
**Marrvel** *                  auto         `⬇ <#marrvel>`_
**CassaNatGenet2017** *        manual       `⬇ <#cassa>`_
**GTEx** *                     auto         `⬇ <#gtex>`_
**BrainSpan** *                auto         `⬇ <#brainspan>`_
**BrainAtlas** *               auto         `⬇ <#brainatlas>`_
**GenCode** *                  auto         `⬇ <#genecode>`_
============================== ============ ===============

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
      :emphasize-lines: 3,6
    
      mutanno download \
              -source_path datasource_directory \
              -source vep \
              -version latest \
              -refversion hg38 \
              -websource mutanno

Make VEP result files and then run preprocess
*********************************************

1. make mock vcf files

   .. code-block::
      
      mutanno vcfmaker \
              -out test.vcf
         

2. run VEP
   

   .. code-block::
      
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


   .. note:: 
   
      `How to install VEP. <https://uswest.ensembl.org/info/docs/tools/vep/script/vep_download.html>`_
   
      1. download `VEP file (v99) from ENSEMBL <ftp://ftp.ensembl.org/pub/release-99/variation/vep/homo_sapiens_vep_99_GRCh38.tar.gz>`_ .
      2. untar and ungzip the downloaded file. 
      3. install plugins (https://uswest.ensembl.org/info/docs/tools/vep/script/vep_plugins.html)


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

For the gene annotation, MutAnno requires several annotation data files from the public resources. You can download and preproces each data sources, and download the preprocessed files from MutAnno storage using the followin comand.

  .. code-block::
     :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source geneannot \
              -version latest \
              -refversion hg38
              -websource mutanno


ENSEMBLgene
^^^^^^^^^^^

* web resource: `ENSEMBL ftp <ftp://ftp.ensembl.org/pub/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.

  .. code-block::
     :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source ensemblgene \
              -version latest \
              -refversion hg38

ENSEMBLgeneGRCh37
^^^^^^^^^^^^^^^^^

* web resource: `ENSEMBL ftp <ftp://ftp.ensembl.org/pub/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.

  .. code-block::
     :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source ensemblgenegrch37 \
              -version latest \
              -refversion hg38

CYTOBAND
^^^^^^^^

* web resource: `UCSC ftp <http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.

  .. code-block::
     :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source cytoband \
              -version latest \
              -refversion hg38

RefSeq
^^^^^^

* web resource: `RefSeqGene data from NCBI ftp <ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.

  .. code-block::
     :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source refseq \
              -version latest \
              -refversion hg38

HGNC
^^^^

* web resource: `HGNC data from EBI ftp <ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


  .. code-block::
     :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source hgnc \
              -version latest \
              -refversion hg38

ClinGen
^^^^^^^

* web resource: `ClinGen curation status site <https://search.clinicalgenome.org/kb/curations/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


  .. code-block::
     :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source clingen \
              -version latest \
              -refversion hg38

ClinGenDisease
^^^^^^^^^^^^^^

* web resource: `ClinGen gene-validity site <https://search.clinicalgenome.org/kb/gene-validity>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


  .. code-block::
     :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source clingendisease \
              -version latest \
              -refversion hg38

ENSEMBLIDxrefTrscriptID
^^^^^^^^^^^^^^^^^^^^^^^

* web resource: `Uniprot ftp <ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


  .. code-block::
     :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source ensemblidxreftrscriptid \
              -version latest \
              -refversion hg38

ENSEMBLIDxref
^^^^^^^^^^^^^

* web resource: `ENSEMBL ftp <ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


  .. code-block::
     :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source ensemblifxref \
              -version latest \
              -refversion hg38

dbNSFP
^^^^^^

* web resource: `dbNSFP ftp <ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


  .. code-block::
     :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source dbnsfp \
              -version latest \
              -refversion hg38

gnomADmetrics
^^^^^^^^^^^^^

* web resource: `gnomAD <ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


  .. code-block::
     :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source gnomadmetrics \
              -version latest \
              -refversion hg38

Marrvel
^^^^^^^

* web resource: `Marrvel <ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.

  .. code-block::
     :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source marrvel \
              -version latest \
              -refversion hg38

CassaNatGenet2017
^^^^^^^^^^^^^^^^^

* web resource: `s_het score from Cassa et al. Nat. Genet. 2017 <https://www.biorxiv.org/highwire/filestream/20869/field_highwire_adjunct_files/0/075523-1.xlsx>`_
* MutAnno doesn't support to download raw source file and preprocess automatically. But it supports to download the preprocessed file from MutAnno storage.

  .. code-block::
     :emphasize-lines: 3,6

      mutanno download \
              -source_path datasource_directory \
              -source ccassanatgenet2017 \
              -version latest \
              -refversion hg38
              -websource mutanno
  


GTEx
^^^^

* web resource: `GTEx dataset <https://gtexportal.org/home/datasets>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


  .. code-block::
     :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source gtex \
              -version latest \
              -refversion hg38

BrainSpan
^^^^^^^^^

* web resource: `BrainSpan ftp <ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


  .. code-block::
     :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source brainspan \
              -version latest \
              -refversion hg38

BrainAtlas
^^^^^^^^^^

* web resource: `BrainAtlas ftp <ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.

  .. code-block::
     :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source brainatlas \
              -version latest \
              -refversion hg38

GenCode
^^^^^^^

* web resource: `GenCode ftp <ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/>`_
* MutAnno supports to 1) download and preprocess automatically, 2) download preprocessed files from MutAnno dropbox, 3) download manually and then run preporcess moduels.


  .. code-block::
     :emphasize-lines: 3
    
      mutanno download \
              -source_path datasource_directory \
              -source gencode \
              -version latest \
              -refversion hg38

