==========
KMOS tools
==========

.. image:: https://img.shields.io/travis/charlottenosam/kmos_tools.svg
        :target: https://travis-ci.org/charlottenosam/kmos_tools


Python tools developed by Charlotte Mason for post-processing and analysing `VLT/KMOS <https://www.eso.org/sci/facilities/paranal/instruments/kmos.html>`_ data from the KLASS VLT Large Program (`196.A-0778 <https://ui.adsabs.harvard.edu/abs/2019Msngr.176...33F/abstract>`_

**Documentation: https://charlottenosam.github.io/kmos_tools.**

Features
--------

* Additional sky subtraction to individual frames
* Measuring shifts between frames from star profiles

Installation
-------------
Navigate to the directory you want to install the code in and clone the `github repository <https://github.com/charlottenosam/kmos_tools>`_:

.. code-block:: bash

    $ git clone git@github.com:charlottenosam/kmos_tools.git
    $ cd kmos_tools
    $ python setup.py install

References and Acknowledgements
-------------------------------

If you find these useful please cite the KLASS survey papers: `Mason et al. (2017) <https://ui.adsabs.harvard.edu/abs/2017ApJ...838...14M/abstract>`_ and `Mason et al. (2019) <https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.3947M/abstract>`_ where I describe some of these methods.

I thank Trevor Mendel, Owen Turner and the ESO helpdesk and KMOS team for many tips. I also benefited from reading the excellent discussion of sky residual cleaning in `Stott et al. (2016) <https://ui.adsabs.harvard.edu/abs/2016MNRAS.457.1888S/abstract>`_ which others may find helpful. I made this package and documentation with a lot of help from `this great tutorial <https://nsls-ii.github.io/scientific-python-cookiecutter/index.html>`_.
