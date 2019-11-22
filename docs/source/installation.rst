============
Installation
============

Navigate to the directory you want to install the code in and clone the `github repository <https://github.com/charlottenosam/kmos_tools>`_:

.. code-block:: bash

    $ git clone git@github.com:charlottenosam/kmos_tools.git
    $ cd kmos_tools
    $ python setup.py install

Requirements
-------------

ESO pipeline
^^^^^^^^^^^^
- ESO Reflex
- `ESO KMOS pipeline <http://www.eso.org/sci/software/pipelines/kmos/kmos-pipe-recipes.html>`_ at least v.1.4

Python packages
^^^^^^^^^^^^^^^
- astropy
- lmfit


You should then be able import in python::
	
	import kmos_tools as kt