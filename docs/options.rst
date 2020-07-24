Optional arguments
==================

.. code::

   mutanno [sub-command] [options]


optional arguments:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit

**sub-commands:**

============== =========================================
``annot``      `run annotation <#run-annotation-annot>`_
``makedata``   `make a single data source file <#make-a-single-data-source-file-makedata>`_
``convert``    `convert <l#convert-file-format-convert>`_
``precal``     pre-calculate
``validate``   validate data format
``preprocess`` quality metrics for VCF
``web``        web view
============== =========================================



Run annotation (``annot``)
--------------------------

  -h, --help            show this help message and exit
  -vcf              
                        input VCF file
  -out              
                        title of output file
  -outtype [vcf, json]
                        output file type
  -ds             
                        data structure json file
  -sourcefile
                        data source file
  -genoinfo
                        add genotype info. in INFO field
  -hgvs                 add hgvs
  -variant_class        add variant class
  -hg19                 add hg19 coordinates
  -chain          
                        chain file for liftover of hg19 coordinates
  -genetable            add gene table
  -blocksize  
                        block size
  -split_multi_allelic_variant
                        split multi-allelic variants
  -clean_tag
                        remove previous annotation information
  -single_source_mode   
                        single source mode
  -load_source_in_memory
                        loading data source in memory
  -sparse               for sparse variant position
  -log
                        log file
  -silence              do not print any log.
  -debug                turn on the debugging mode



Make a single data source file (``makedata``)
---------------------------------------------


Convert file format (``convert``)
---------------------------------


