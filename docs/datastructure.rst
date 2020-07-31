Data Structure File
===================

If you download data source files with ``mutanno download`` command, **MutAnno** generates data structure file (``.json``) automatically.

If you want to add or remove some fields, you can edit this file.

Basic strcture
--------------

The basic heirachical structure of json file is:

.. code::

   {
      "name": "MUTANNO",
      ...
      "source": [
        {
        	"name": "SOURCE1",
        	...
        	"fields": [
        		{
        			"name": "FIELD1",
        			...
        		}
        	]
        }
      ]
   }



Root attributes
---------------

* ``name`` : root name (default: "MUTANNO")
* ``desc`` : description 
* ``version`` : MutAnno version
* ``version_date`` : MutAnno version date
* ``level`` : annotation level (default: "variant")
* ``sourcefile_path`` : source file path(directory)
* ``sourcefile`` : single source file path

Source attributes
-----------------

* ``name`` : source name (default: "MUTANNO")
* ``desc`` : description 
* ``version`` : source version
* ``version_date`` : source version date
* ``subembedded`` : if this source subembedded fields. 
* ``function`` : TBA
* ``params`` : TBA
* ``is_available`` : is avilable soruce? (`true`, `false`; default=`true`)
* ``sourcefile`` : this source file path
* ``format`` : source file format (`tsi`, `bed`, `vcf`; default=`tsi`)

Field attributes
----------------

* ``name`` : field name (default: "MUTANNO")
* ``desc`` : description 
* ``type`` : field variable type (`string`, `boolean`, `number`) 
* ``is_list`` : is list type field (`true`, `false`; default=`false`)
* ``delimiter`` : if the field type is `list` (which means ``is_list: true``), the delimiter can be assigned.
* ``function`` : TBA
* ``params`` : TBA
* ``is_available`` : is avilable field? (`true`, `false`; default=`true`)