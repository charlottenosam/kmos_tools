# -*- coding: utf-8 -*-
"""Class to hold individual KMOS exposures for post-processing

Takes ``*SCI_RECONSTRUCTED.fits`` filenames as input

Usage:
    Create a instance of the exposure class for each exposure file::

        from kmos_tools import exposure

        fname = 'KMOS_SCI_RECONSTRUCTED.fits'
        exp = exposure.Exposure(fname)

"""

__author__ = ["Charlotte Mason (UCLA)", "Antonello Calabro' (OAR)"]

import astropy.io.fits as fits


class Exposure(object):
    """A single KMOS exposure object

    Takes SCI_RECONSTRUCTED filename as input, to be used for other operations.

    Attributes:
        filename (str): filepath to SCI_RECONSTRUCTED fits file
        hdulist (astropy sequence of HDU objects or single HDU): HDUlist of IFUs in the exposure
        hdr (astropy header): FITS header
        filter (str): KMOS filter used
        frame_time (str): date of observations


    """
    def __init__(self, reconstructed_fits_path, vb=False):
        """init the Exposure class using SCI_RECONSTRUCTED files

        Args:
            reconstructed_fits_path (str): File path of SCI_RECONSTRUCTED individual frame.
            vb (bool): verbose?


        """

        self.filename = reconstructed_fits_path
        self.hdulist  = fits.open(self.filename)
        self.hdr      = self.hdulist[0].header
        self.filter   = self.hdr['HIERARCH ESO INS FILT1 ID']

        self.frame_time = self.hdr['DATE-OBS']

        self.vb = vb

        if self.vb:
            print('REDUCTION: Inspecting %s' % self.filename)

        # Work out whether we are in science (A) or sky (B) mode
        count = 0
        for ifu in range(1, 25):
            try:
                ext = self.hdulist['IFU.%i.DATA' % ifu]
                if ext.data is not None:
                    count = count + 1
            except:
                if self.vb:
                    print('WARNING: no IFU %i' % ifu)

        if count > 5:
            self.mode = 'A'
        else:
            self.mode = 'B'

        self.star_ifu = None
        self.filename_fluxfix = self.filename_fluxfix = self.filename.replace('.fits', '_FLUXFIX.fits')
        self.filename_skycorr = self.filename.replace('.fits', '_SKYSUB.fits')
        self.starfile = None

        return
