Make Single Source File (``makedata``)
======================================

After downloading data source files, basically you don't have to make the single source file. But you want to manage the source files more efficiently, you can make a single source file from serveral source files.


.. code::

   mutanno makedata \
           -ds data_structure.json \
           -out single_datasource_file.tsi \
           -vartype SNV \
           -blocksize 1000


Input file
----------

To make single source file, **MutAnno** requires a data structure file (``-ds``) and data source files that the structure file contains.

Data structure file
^^^^^^^^^^^^^^^^^^^

The data structure file must include data source file path. 


Data source file
^^^^^^^^^^^^^^^^

You can download various data source files. Please see `here <download.html>`_. Each source file path should be inserted in the structure file.


.. code:: json

   {
     "name": "MUTANNO",
     "source": [
       {
         "name": "DBSNP",
         "version": "20190722",
         "version_date": "07/22/2019",
         "sourcefile": "ANNOTATION_SOURCE/GRCh38_latest_dbSNP_all.convchrom.onlyid.vcf.gz",
         "format": "tsi"
       }
      ]
    }


Output file
-----------

With ``-out`` option, you can assign the output file name and path for the single source file.



Annotation type
---------------

By default, ``-vartype`` is ``SNV`` among ``SNV``, ``GENE``, and ``CODING_GENE_MAIN_CHROM``. 

SNV
^^^

GENE
^^^^

To make gene anntoation file, we can use


CODING_GENE_MAIN_CHROM
^^^^^^^^^^^^^^^^^^^^^^





Only for specific genomic region
--------------------------------

If you want to make the single annotation file for specific genomic region, you can use ``-region`` option.

.. code::

   mutanno makedata \
           -ds data_structure.json \
           -out single_datasource_file.tsi \
           -vartype SNV \
           -blocksize 1000
           -region chr1:1000000-2000000


