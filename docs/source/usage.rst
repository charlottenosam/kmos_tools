====================
KMOS Reduction Steps
====================

These are my steps for reducing KMOS data

1. Download data from ESO
=========================
Unzip the data::

	$ gzip -d *.Z

If you want to inspect the final data you can use::

	$ dfits –x 0 KMOS*fits | fitsort tpl.id det.seq1.dit ins.filt1.id ocs.rot.naangle

.. note:: May need to correct bias from detector readout (causes striping in raw images), with code from Trevor Mendel (but I think this is fixed in recent versions of the ESO pipeline)


2. Run EsoReflex
================
.. code-block:: bash

	esoreflex &

I use the following settings:

- Remove the COMBINE step so that all OBs can be combined later
- SCI_RED ==> sky_tweak = true
- Set UseSkyFlats = no (to use lamp flats instead - taken at more rotation angles, better according to ESO)
- Animate at runtime 100ms

3. Sky subtraction correction
=============================

Do additional sky subtraction correction to individual exposures

.. note:: ESO pipeline v<1.4 may need to do flux calibration correction, to divide by DIT length. Can use ``kmos_tools.pipeline_fixes.flux_fix()``.

.. code-block:: python

	import kmos_tools as kt

	# Find your individual exposures
	exp_list = kt.find_exposures(prefix='*', dir='*/')

	for exp in exp_list:
		frame = kt.Exposure(exp)

	    # Make 1D sky residual spectra for each exposure (1 per detector)
	    kt.make_sky_residual_spectra(frame)

	    # Do additional sky subtraction and save
	    kt.subtract_sky_residual_spectra(frame)

.. note:: I also had some success with `ZAP <https://zap.readthedocs.io/en/2.1/>`_. however it had problems at the edges of the cubes.


4. Rotate frames which don't have North up (optional)
=====================================================

Sometimes the instrument is rotated between OBs, need to rotate the IFUs so they can be combined

.. code-block:: python

	import kmos_tools as kt

	# Find your individual exposures
	exp_list = kt.find_exposures(prefix='*', dir='*/')

	for exp in exp_list:
		frame = kt.Exposure(exp)

	    kt.rotation_fix(frame)


5. Measure star PSFs and calculate user shifts
==============================================

Use observed stars in each exposure to estimate the PSF and then calculate a list of accurate dither offsets between exposures (as KMOS doesn't always go where you set the dithers...)

5.1. Find positions of stars
-----------------------------

Finds stars in exposures, measures their PSF FWHM and positions and creates a list of the 'good' frames (below a PSF cut) and a table of the star positions.

.. code-block:: python

	# List of frames to check (should be all the frames you want to combine)
	frame_list = ['frame1.fits', 'frame2.fits']

	fname_stars   = 'stars.txt'   # File to save star parameters to
	fname_combine = 'combine.sof' # .sof file to create

	kt.star_positions_batch(frame_list, psf_cut=0.8, starparams_filename=fname_stars, combinefiles_filename=fname_combine)

5.2. Calculate shifts between frames
------------------------------------

Using the ``fname_stars`` file, calculate the shifts between each frame.

.. code-block:: python

	fname_usershifts = 'usershifts.txt' # file to save usershifts to

	kt.make_user_shifts_file(starparams_filename=fname_stars, usershifts_filename=fname_usershifts)

.. note:: If using ``user`` shifts to combine KMOS frames, the order of frames in ``.sof`` **must match** the order of the offsets in ``'usershifts.txt'``. This should be done automatically in these scripts, so don't reorder the frames in the ``.sof`` file.

6. Combine exposures
====================

I use the default (sigma clipping) combination, using my calculated shifts between frames

.. note:: In v1.4 combine.sof file needs SCI_RECONSTRUCTED after filename

.. code-block:: bash

	combinedir=combineOBs
	mkdir $combinedir
	cp usershifts.txt combine.sof $combinedir
	cd $combinedir
	esorex kmos_combine ==method='user' ==filename='usershifts.txt' ==edge_nan combine.sof
