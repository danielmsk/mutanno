Download Data Source (``download``)
===================================

To annotate variants, **MutAnno** need data source files and can download them with ``download`` sub-command.

.. code::

   mutanno download \
           -source_path datasource_directory \
           -source all \
           -version lastest \
           -refversion hg38


* ``-source_path`` : data source path (default: `mutanno_source`)
* ``-source`` : source name (default: `all`)
* ``-version`` : source version (default: `lastest`)
* ``-refversion`` : reference version [`hg19`, `hg38`] (default: `hg38`)

.. note::

   If ``-source`` is ``all``, all data sources can be downloaded and preprocessed.

Variant database
----------------

ClinVar
^^^^^^^

.. code::

   mutanno download \
           -source_path datasource_directory \
           -source CLINVAR \
           -version lastest \
           -refversion hg38

MutAnno can download CLINVAR data source files (.vcf and .xml) from NCBI FTP site and run pre-processing automatically. 
Finally it generates mti files as the following directory structure.

.. code::

   datasource_directory/
   |- downloaded/CLINVAR/[hg38/hg19]/[version]/
   |  |- clinvar_OOOOOO.vcf.gz
   |  |- ClinVarFullRelease_OOOOOO.xml.gz
   |- CLINVAR_hg38_[version]_submission.mti.gz
   |- CLINVAR_hg38_[version]_variant.mti.gz


Internal working flow
---------------------

1. :guilabel:`/src/mutanno/__init__.py` : dispatch job

   a. ``cli()``
   b. ``dispatch_job()`` -> call ``Downloader.run()``

2. :guilabel:`/src/mutanno/download.py` : download data source file

   a. ``Downloader.run()``
   b. ``Downloader.load_sourcelist_file()`` : load :guilabel:`/src/mutanno/data/sourcelist.json`
   c. ``Downloader.download_and_preprocess_sourcefile()``

      i. ``Downloader.download_source_file()`` : download source files from public sites.
      ii. ``Downloader.preprocess_source_file()`` : call pre-processing function (:guilabel:`/src/mutanno/preprocess/[designated_function]`)

         * The designated_function is from :guilabel:`/src/mutanno/data/sourcelist.json`.

            .. code-block:: json
               :emphasize-lines: 8

               {
               "title":"MUTANNO_SOURCE_LIST",
               "sources":{
                  "hg38":[
                     {
                     "name":"CLINVAR",
                     "latest": "20200329",
                     "preprocess_function": "run_preproc_clinvar",
                     "versions":[

         * This designated_function is defined in :guilabel:`/src/mutanno/preprocess/__init__.py`.