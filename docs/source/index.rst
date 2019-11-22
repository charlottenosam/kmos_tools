.. Packaging Scientific Python documentation master file, created by
   sphinx-quickstart on Thu Jun 28 12:35:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

KMOS tools Documentation
========================

These are some key python routines used to post-process `VLT/KMOS <https://www.eso.org/sci/facilities/paranal/instruments/kmos.html>`_ data in the KLASS VLT Large Program (`196.A-0778 <https://ui.adsabs.harvard.edu/abs/2019Msngr.176...33F/abstract>`_).

This not complete, tested on other systems, or very well-supported... If you have problems I encourage you to fix them and submit a `pull request <https://github.com/charlottenosam/kmos_tools/pulls>`_ :)!


References and Acknowledgements
-------------------------------

If you find these useful please cite the KLASS survey papers: `Mason et al. (2017) <https://ui.adsabs.harvard.edu/abs/2017ApJ...838...14M/abstract>`_ and `Mason et al. (2019) <https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.3947M/abstract>`_ where I describe some of these methods.

I thank Trevor Mendel, Owen Turner and the ESO helpdesk and KMOS team for many tips. I also benefited from reading the excellent discussion of sky residual cleaning in `Stott et al. (2016) <https://ui.adsabs.harvard.edu/abs/2016MNRAS.457.1888S/abstract>`_ which others may find helpful. I made this package and documentation with a lot of help from `this great tutorial <https://nsls-ii.github.io/scientific-python-cookiecutter/index.html>`_.


Contents
---------

.. toctree::
   :maxdepth: 2

   installation
   usage
   api
   release-history
   