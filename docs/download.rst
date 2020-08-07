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
           -source clinvar \
           -version lastest \
           -refversion hg38
