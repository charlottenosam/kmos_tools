# kmos_tools/sky_clean.py
# -*- coding: utf-8 -*-
"""Module to clean up sky subtraction residuals in KMOS data

Works on individual exposures which must then be combined via ``kmos_combine``.


"""

import os
import numpy as np
import astropy.io.fits as fits
import itertools as it
import datetime
import matplotlib.pyplot as plt
from lmfit import Parameters, minimize


__all__ = ["make_sky_residual_spectra", "subtract_sky_residual_spectra", "sky_residual"]


def make_sky_residual_spectra(exposure, clobber=True, plot=True):
    """Calculate all the sky residual corrections for one DIT ~ a la Trevor Mendel

    Generate median 1D sky spectra for each detector,
    only from S1 and S3 cubes, and save them to use for subtraction.

    Args:
        exposure (object): exposure object
        clobber (bool): Make a new sky spectrum if it already exists
        plot (bool): make and save a plot of the median sky spectra

    Returns:
        Saves FITS file with 1D sky spectra for each detector


    """

    print('Correcting sky residuals in', exposure.filename)

    # Create stacks of 'empty' cubes
    exposure.filename_skyspec = exposure.filename.replace('.fits', '_SKYSPEC.fits')
    if os.path.exists(exposure.filename_skyspec) and clobber is False:
        print(exposure.filename_skyspec, 'sky spectra already exist, moving on...')
    else:
        detector1, detector2, detector3 = [], [], []
        for ifu in range(1, 25):

            ext = exposure.hdulist['IFU.' + str(ifu) + '.DATA']

            if len(ext.shape) > 0:

                ifu_comment = 'HIERARCH ESO OCS ARM' + str(ifu) + ' NAME'
                ifu_header  = ext.header

                # Use only frames with S1 or S3 targets for sky subtraction
                if 'S1' in exposure.hdr[ifu_comment] or 'S3' in exposure.hdr[ifu_comment]:

                    ifu_cube = ext.data

                    if 1 <= ifu <= 8:
                        detector1.append(ifu_cube)
                    elif 9 <= ifu <= 16:
                        detector2.append(ifu_cube)
                    else:
                        detector3.append(ifu_cube)

        # Stack all spectra
        detector1, detector2, detector3 = np.array(detector1), np.array(detector2), np.array(detector3)
        if detector1.size > 0 and detector2.size > 0 and detector3.size > 0:
            detector_all    = np.concatenate((detector1, detector2, detector3), axis=0)
        elif detector1.size > 0 and detector2.size > 0:
            detector_all = np.concatenate((detector1, detector2), axis=0)
        elif detector1.size > 0 and detector3.size > 0:
            detector_all    = np.concatenate((detector1, detector3), axis=0)
        elif detector2.size > 0 and detector3.size > 0:
            detector_all    = np.concatenate((detector2, detector3), axis=0)
        else:
            if detector1.size > 0:
                detector_all = detector1
            elif detector2.size > 0:
                detector_all = detector2
            elif detector3.size > 0:
                detector_all = detector3

        # Generate median of 'empty' stacks to use as sky residual spectra for each detector
        skyspec_1D_all = np.nanmedian(detector_all, axis=(0, 2, 3))
        skyspec_1D = {}
        detectors = [detector1, detector2, detector3]
        for i in range(len(detectors)):
            if detectors[i].shape[0] > 1:
                skyspec_1D[i] = np.nanmedian(detectors[i], axis=(0, 2, 3))
            else:
                skyspec_1D[i] = None

        if plot:
            plt.figure(figsize=(10, 5))

            plt.plot(skyspec_1D_all, lw=1, alpha=0.8, label='All detectors (%i IFUs)' % detector_all.shape[0], zorder=10)

            for i in range(len(skyspec_1D)):
                if skyspec_1D[i] is not None:
                    plt.plot(skyspec_1D[i], lw=1, alpha=0.8, label='Detector %i (%i IFUs)' % (i, detectors[i].shape[0]))

            ymin, ymax = np.nanpercentile(skyspec_1D_all, 1), np.nanpercentile(skyspec_1D_all, 99)
            if ymin > 0.: ymin = -1.e-18

            plt.ylim(ymin, ymax)
            plt.xlabel('Wavelength [pix]')
            plt.ylabel('Flux')
            plt.title('Sky Subtraction Residuals - from median S1 and S3')
            plt.legend()
            plt.tight_layout()

            plt.savefig(exposure.filename_skyspec.replace('.fits', '.pdf'))

        # Save spectra to fits file
        # Headers
        prihdr = exposure.hdr.copy()
        prihdr.add_comment('Sky Spectrum on each detector, and median sky spectrum')
        hdr1D = fits.Header()
        hdr1D['SIMPLE'] = 'T'
        hdr1D['BITPIX'] = -32
        hdr1D['NAXIS'] = 1
        hdr1D['NAXIS1'] = 2048
        hdr1D['PCOUNT'] = 0
        hdr1D['GCOUNT'] = 1
        hdr1D['CUNIT1'] = ifu_header['CUNIT3']
        hdr1D['CRPIX1'] = ifu_header['CRPIX3']
        hdr1D['CRVAL1'] = ifu_header['CRVAL3']
        hdr1D['CDELT1'] = ifu_header['CDELT3']
        hdr1D['BUNIT'] = 'cgs'

        hdr1D_1, hdr1D_2, hdr1D_3 = hdr1D.copy(), hdr1D.copy(), hdr1D.copy()
        hdr1D['EXTNAME']   = 'ALL'
        hdr1D_1['EXTNAME'] = 'DETECTOR1'
        hdr1D_2['EXTNAME'] = 'DETECTOR2'
        hdr1D_3['EXTNAME'] = 'DETECTOR3'

        # Extensions
        hdu     = fits.PrimaryHDU(header=prihdr)
        hdu_all = fits.ImageHDU(skyspec_1D_all, header=hdr1D)
        hdu_1   = fits.ImageHDU(skyspec_1D[0], header=hdr1D_1)
        hdu_2   = fits.ImageHDU(skyspec_1D[1], header=hdr1D_2)
        hdu_3   = fits.ImageHDU(skyspec_1D[2], header=hdr1D_3)

        # Create hdu list and write
        hdulist = fits.HDUList([hdu, hdu_all, hdu_1, hdu_2, hdu_3])
        hdulist.writeto(exposure.filename_skyspec, clobber=True)
        print('Saved fits file to ', exposure.filename_skyspec)

    return


def subtract_sky_residual_spectra(exposure, clobber=True, plot=True):
    r"""Subtract residual sky spectrum from each IFU using relevant detector sky spectrum

    Rescales 1D sky spectrum :math:`s_\lambda` in each spatial pixel by :math:`A_{x, y}` \
    for each IFU such that sky subtracted spectrum:

    .. math:: f_{x, y, \lambda}^\mathrm{skysub} = f_{x, y, \lambda} - A_{x, y}s_\lambda

    Where :math:`A_{x, y}` is obtained by minimising:

    .. math:: \left( \frac{f_{x, y, \lambda}^\mathrm{skysub}}{\sigma_\lambda} \right)^2

    Args:
        exposure (object): exposure object
        clobber (bool): Make a new sky spectrum if it already exists
        plot (bool): make and save a plot of the median sky spectra

    Returns:
        Saves FITS file with 1D sky spectra for each detector with filename ``*_SKYSUB.fits``


    """

    # Do sky subtraction
    if os.path.exists(exposure.filename_skycorr) and clobber is False:

        print(exposure.filename_skycorr, 'already exists, moving on...')

    else:

        skyspec_all = fits.open(exposure.filename_skyspec)

        for ifu in range(1, 25):

            ext = exposure.hdulist['IFU.' + str(ifu) + '.DATA']

            if len(ext.shape) > 0:

                # Finding detector
                if 1 <= ifu <= 8:
                    detector = 1
                elif 9 <= ifu <= 16:
                    detector = 2
                elif 17 <= ifu <= 24:
                    detector = 3

                skyspec  = skyspec_all[detector+1].data

                # Do sky subtraction
                if skyspec is None:
                    print('No sky spectrum for detector %i, so using the median' % detector)
                    skyspec = skyspec_all[1].data

                # Estimate 1D error from std of flux
                dat3D = ext.data.copy()
                err1D = np.nanstd(dat3D, axis=(1, 2))/np.sqrt(dat3D.shape[1] * dat3D.shape[2])

                # Get rescaling # TODO better minimise
                skyscale2D   = np.zeros_like(dat3D[0])
                chisquared2D = np.zeros_like(dat3D[0])
                for i, j in it.product(np.arange(dat3D.shape[1]), np.arange(dat3D.shape[2])):
                    params = Parameters()
                    params.add('scale', value=1.)
                    out = minimize(exposure.sky_residual, params, args=(skyspec, dat3D[:, i, j], err1D))
                    skyscale2D[i, j]   = out.params['scale'].value
                    chisquared2D[i, j] = out.chisqr

                # Sky subtract
                data_corr = dat3D - skyscale2D*skyspec[:, None, None]

                # Subtract median flux if not S2
                ifu_comment = 'HIERARCH ESO OCS ARM' + str(ifu) + ' NAME'
                if 'S2' not in exposure.hdr[ifu_comment]:
                    data_corr -= np.nanmedian(data_corr)

                ext.data  = data_corr
                now       = datetime.datetime.now()
                ext.header['SKY RESIDUALS CORRECTED'] = str(now)

        exposure.hdulist.writeto(exposure.filename_skycorr, clobber=clobber)
        print('Saved Sky Fix to', exposure.filename_skycorr)

        return


def sky_residual(params, sky, data, err):
    r"""Residual for 1D sky-scaling

    Used to find optimal rescaling of 1D spectrum, A,
    for each spatial pixel

    .. math::

        R_i = \frac{(f_i - A s_i)}{\sigma_i} \\
        \chi^2 = \sum_i R_i^2


    Args:
        params (params): parameters
        sky (ndarray): 1D sky spectrum
        data (ndarray): 3D data cube
        err (ndarray): 1D error spectrum

    Returns:
        residual


    """
    scale = params['scale']
    return np.nan_to_num((data - scale * sky)/err)
