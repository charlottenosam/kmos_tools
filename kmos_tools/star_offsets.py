# kmos_tools/star_offsets.py
# -*- coding: utf-8 -*-
"""Module to find stars in an exposure, and fit their spatial positions to use in ``kmos_combine --user_shifts``


    Usage::

        # List of frames to check (should be all the frames you want to combine)
        frame_list = ['frame1.fits', 'frame2.fits']

        fname_stars   = 'stars.txt'   # File to save star parameters to
        fname_combine = 'combine.sof' # .sof file to create

        kt.star_positions_batch(frame_list, psf_cut=0.8, starparams_filename=fname_stars, combinefiles_filename=fname_combine)

        fname_usershifts = 'usershifts.txt' # file to save usershifts to

        kt.make_user_shifts_file(starparams_filename=fname_stars, usershifts_filename=fname_usershifts)

"""

import os
import numpy as np
import subprocess
import astropy.io.fits as fits
import datetime
import kmos_tools as kt


__all__ = ["make_user_shifts_file", "star_positions_batch", "find_star_ifu", "star_fit_profile", "star_psf"]


def make_user_shifts_file(starparams_filename, usershifts_filename=None):
    """Create file of user shifts to use with ``kmos_combine --method="user"``

    Note:
        To run ``kmos_combine`` you *must* use the frames in the order they are used here
        i.e. so the offset matches the right frame

    Args:
        starparams_filename (str): filepath of table with the star positions
        usershifts_filename (str, optional): filename to save usershift table to


    """

    # Get info from star table
    star_table = np.genfromtxt(starparams_filename, dtype=None, names=True, skip_header=1)
    n_frames   = len(star_table)

    # Calculate shifts and save to file
    print('Calculating shifts')
    shiftx, shifty = [0]*(n_frames-1), [0]*(n_frames-1)
    for i in range(1, n_frames):
        shiftx[i-1] = star_table['XCEN_pix'][i] - star_table['XCEN_pix'][0]
        shifty[i-1] = star_table['YCEN_pix'][0] - star_table['YCEN_pix'][i]

    if usershifts_filename is None:
        usershifts_filename = starparams_filename.replace('.txt', '_usershifts.txt')

    np.savetxt(usershifts_filename, np.array([shiftx, shifty]).T, fmt='%f', delimiter='    ')
    print(' - Saved shifts file to', usershifts_filename)

    return


def star_positions_batch(frame_list, psf_cut=0.8, edge_x=2., edge_y=2., star_ifu=None,
                         starparams_filename=None, combinefiles_filename=None):
    """Given a list of exposure file names, will find stars and measure PSFs.

    Given a list of exposure file names, will find stars and measure PSFs.
    Then creates a list of frames with PSF FWHM below a given value, saves
    the list of star positions (to calculate user shifts) and creates a .sof
    file of the 'good' frames.

    Args:
        frame_list (list): List of file names of individual ``SCI_RECONSTRUCTED.fits`` files
        psf_cut (float): max PSF to use for combining (in arcsec) [default = 0.8 arcsec]
        edge_x (float): x_star > edge_x to include (i.e. not on the edge) [default = 2 pixels]
        edge_y (float): y_star > edge_y to include (i.e. not on the edge) [default = 2 pixels]
        star_if (int): IFU star is on

    Yields:
        starparams_filename (str): table with parameters of good stars
        combinefiles_filename (str): .sof file with list of good frames

    """
    combine_files = []
    star_table, star_table_bad = [], []

    if starparams_filename is None:
        starparams_filename = frame_list[0].split('_products')[0]+'_products/star_psf_'+str(datetime.date.today())+'.txt'
    if combinefiles_filename is None:
        combinefiles_filename = frame_list[0].split('_products')[0]+'_products/combine_'+str(datetime.date.today())+'.sof'

    for ff, frame in enumerate(frame_list):
        sci_reconstructed = kt.Exposure(frame)
        sci_reconstructed.star_ifu = star_ifu

        # Look for a star in the frame (based on having `star` in the target name) and fit a gaussian profile to it
        try:
            psf_center_x, psf_center_y, psf_fwhm, psf_ba, psf_pa, invert_comment = star_psf(sci_reconstructed,clobber=True)
            print(psf_center_x)
        except:
            print('No star in %s' % frame)
            psf_center_x, psf_center_y, psf_fwhm, psf_ba, psf_pa, invert_comment = np.nan, np.nan, np.nan, np.nan, np.nan, '# no star'

        if psf_fwhm < psf_cut and sci_reconstructed.mode == 'A' and psf_center_x > edge_x and psf_center_y > edge_y:
            combine_files.append([frame, 'SCI_RECONSTRUCTED'])
            star_table.append([frame, sci_reconstructed.frame_time, psf_center_x, psf_center_y, psf_fwhm, psf_ba, psf_pa, invert_comment])
        elif psf_center_x == np.nan:
            star_table_bad.append([frame, sci_reconstructed.frame_time, psf_center_x, psf_center_y, psf_fwhm, psf_ba, psf_pa, invert_comment])
        else:
            print('WARNING: Something wrong with PSF in %s' % frame)
            star_table_bad.append([frame, sci_reconstructed.frame_time, psf_center_x, psf_center_y, psf_fwhm, psf_ba, psf_pa, invert_comment])

        del sci_reconstructed

    # Save star parameter output
    star_table = np.array(star_table)

    med_FWHM   = [float(x) for x in star_table[:, 4]]
    med_FWHM   = np.nanmedian([x for x in med_FWHM if x > 0])
    med_ba     = [float(x) for x in star_table[:, 5]]
    med_ba     = np.nanmedian([x for x in med_ba if x > 0])
    med_pa     = [float(x) for x in star_table[:, 6]]
    med_pa     = np.nanmedian([x for x in med_pa if x > 0])

    np.savetxt(starparams_filename, star_table, fmt='%s', delimiter='    ', header='Star PSFs. Median FWHM = %.3f Median BA = %.3f Median PA = %.3f \nFrame    OB_time     XCEN_pix     YCEN_pix    FWHM_arcsec    BA    PA_deg' % (med_FWHM, med_ba, med_pa), comments='# ')
    np.savetxt(starparams_filename.replace('.txt', '_bad.txt'), star_table_bad, fmt='%s', delimiter='    ', header='Bad Star PSFs. \nFrame    OB_time     XCEN_pix     YCEN_pix    FWHM_arcsec    BA    PA_deg', comments='# ')
    print(' - Saved PSF info file to %s ' % starparams_filename)

    # Save .sof file with all the successful frames
    combine_files = np.array(combine_files)
    np.savetxt(combinefiles_filename, combine_files, fmt='%s')
    print(' - Saved combine.sof file to %s ' % combinefiles_filename)

    return


def find_star_ifu(exposure):
    """Find the IFUs in the SCI_RECONSTRUCTED cubes containing a star

    Args:
        exposure (object): exposure object

    Returns:
        exposure.star_ifu (int): IFU containing star


    """
    if exposure.star_ifu is None:

        for ifu in range(1, 25):

            ifu_comment = 'HIERARCH ESO OCS ARM%i NAME' % ifu

            # Look for star in name of IFU targe
            if 'star' not in exposure.hdr[ifu_comment].lower():
                continue
            else:
                print('Star found in IFU %i' % ifu)
                exposure.star_ifu = ifu

    return exposure.star_ifu


def star_fit_profile(exposure):
    """Fit gaussian profiles to star to find PSFs

    Using ESO command line pipeline tools (``esorex``) extract the star from the exposure,
    fit a Gaussian profile to it, and save the fitted star to a new fits file ``star_file_name``.


    Args:
        exposure (object): exposure object

    Returns:
        star_file_name (str): filepath of FITS file with star with PSF fitted to it
        invert (bool): was the star flux weird and inverted?


    """
    #print('Y')
    # Find IFU star is on
    if exposure.star_ifu is None:
        exposure.star_ifu = find_star_ifu(exposure)

    # Make a copy of the IFU with the star
    kmo_copy = 'esorex kmo_copy -x=1 -y=1 -z=1 -xsize=14 -ysize=14 -zsize=2048 -ifu=%s %s' % (str(exposure.star_ifu), exposure.filename)
    subprocess.call(kmo_copy,shell=True)
    status, copyfile = subprocess.getstatusoutput("find . -maxdepth 1 -iname copy.fits")

    # Collapse the IFU to make an image
    kmo_make_image = 'esorex kmo_make_image %s' % copyfile
    subprocess.call(kmo_make_image,shell=True)
    status, makeimagefile = subprocess.getstatusoutput("find . -maxdepth 1 -iname make_image.fits")
    # Check the image, if weird, invert
    image_hdu = fits.open(makeimagefile)
    image     = image_hdu[1].data
    test_star = np.nansum(image[3:-3, 3:-3] - np.nanmedian(image[3:-3, 3:-3]))
    invert = False
    if test_star < 0.:
        print('WARNING: Weird star in %s, multiplying image by -1' % exposure.filename)
        image_hdu[1].data = -1. * image
        image_hdu.writeto(makeimagefile, clobber=True)
        invert = True

    # Rename star file
    exposure.star_image_file = exposure.filename.strip('.fits') + '_star_image.fits'
    subprocess.call('cp %s %s' % (makeimagefile, exposure.star_image_file),shell=True)

    # Fit profile to star image
    kmo_fit_profile = 'esorex kmo_fit_profile %s' % exposure.star_image_file
    subprocess.call(kmo_fit_profile,shell=True)
    status, fitprofilefile = subprocess.getstatusoutput("find . -maxdepth 1 -iname fit_profile.fits")

    # Tidy up
    star_file_name     = exposure.filename.strip('.fits')+'_star_psf.fits'
    rename_fit_profile = 'mv %s %s' % (fitprofilefile, star_file_name)
    delete_temps       = 'rm %s %s' % (copyfile, makeimagefile)

    subprocess.call(rename_fit_profile,shell=True)
    subprocess.call(delete_temps,shell=True)

    print('Saved star psf profile to '+star_file_name)

    # Filename to save star IFU to
    exposure.starfile = exposure.filename.strip('.fits')+'_star_psf.fits'

    return star_file_name, invert


def star_psf(exposure, clobber=True, vb=False):
    """Get the PSF profiles of the star for calculating user shifts

    Read the file ``star_file_name`` and return the main parameters


    Args:
        exposure (object): exposure object

    Returns:
        psf_center_x (float): PSF centroid x in pixels
        psf_center_y (float): PSF centroid y in pixels
        psf_fwhm (float): PSF FWHM in arcsec
        psf_ba (float): PSF BA # TODO I don't remember what this is!!
        psf_pa (float): PSF position angle in degrees
        exposure.invert (bool): is the flux weird?


    """

    exposure.invert = False
    #print('0')
    #if clobber is True and ~os.path.exists(exposure.starfile): - problem with this line
    #print('1')
    exposure.starfile, exposure.invert = star_fit_profile(exposure) #changed this line - wasn't calling function
    #print(exposure.starfile,exposure.invert)

    if exposure.invert:
        invert_comment = '# weird star, inverted flux'
    else:
        invert_comment = ''

    star_hdulist  = fits.open(exposure.starfile)
    star_hdr      = star_hdulist[1].header
    # Load fit parameters
    psf_center_x = star_hdr['HIERARCH ESO PRO FIT CENTROID X']
    psf_center_y = star_hdr['HIERARCH ESO PRO FIT CENTROID Y']
    psf_r_x      = star_hdr['HIERARCH ESO PRO FIT RADIUS X']
    psf_r_y      = star_hdr['HIERARCH ESO PRO FIT RADIUS Y']

    # Get FWHM, BA, and PA
    psf_fwhm = 2.3548 * 0.5 * (psf_r_x + psf_r_y) * 0.2  # fwhm in arcsec
    psf_ba   = psf_r_y / psf_r_x
    psf_pa   = star_hdr['HIERARCH ESO PRO FIT ROT']
    if psf_pa > 0.:
        psf_pa = psf_pa % 360.

    if vb:
        print('XPIX=%.2f, YPIX=%.2f, FWHM=%.3f arcsec, BA=%.3f, PA=%.3f' % (psf_center_x, psf_center_y, psf_fwhm, psf_ba, psf_pa))

    return psf_center_x, psf_center_y, psf_fwhm, psf_ba, psf_pa, invert_comment
