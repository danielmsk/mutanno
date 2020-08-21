Installation
============


Prerequisites
-------------

* python 3.4+
* `pysam (0.11.2.2+) <https://pypi.org/project/pysam/>`_
* `pyfaidx (0.5.3.1+) <https://pypi.org/project/pyfaidx/>`_
* `pytabix (0.0.2+) <https://pypi.org/project/pytabix/>`_
* `pyliftover (0.4+) <https://pypi.org/project/pyliftover/>`_
* `xmltodict (0.12.0+) <https://pypi.org/project/xmltodict/>`_
* `pyvcf (0.6.8+) <https://pypi.org/project/PyVCF/>`_ 
* `tqdm (4.48.2+) <https://pypi.org/project/tqdm/>`_

.. note::
	If you are using pypi (``pip``) for **MutAnno** installation, you don't need to install prerequisite python libraries because pypi can install them automatically. 

Install with pypi
-----------------

To install **MutAnno**, open a shell and run the following command:

.. code::

    $ pip install mutanno


Install with code from github
-----------------------------

Alternatively, use ``git clone`` followed by ``setup.py``

.. code::

    $ git clone https://github.com/dbmi-bgm/mutanno
    $ cd mutanno
    $ python setup.py install

