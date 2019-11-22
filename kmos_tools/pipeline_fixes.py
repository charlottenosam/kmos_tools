# -*- coding: utf-8 -*-
"""Quick fixes for KMOS pipeline problems. Methods are:

    Methods:
        rotation_fix: Rotates frames so that they are all orientated in the same direction for combining
        calib_fix: Fix flux calibration if it wasn't done by the pipeline (i.e. no standard observed on on detector)
        flux_fix: Fixes an old (now fixed??) pipeline problem where flux needs to be divided by DIT length in seconds


"""
import os
import subprocess
import astropy.io.fits as fits
import datetime


__all__ = ["rotation_fix", "flux_fix"]


def rotation_fix(exposure, updateheader=False, keepsize=False):
    """Sometimes the instrument was rotated, this rotates all the frames to the same angle for combining

    For some OBs the field is rotated so North is no longer up. This is when ocs.rot.offangle != 0

    Need to rotate back to North = up so we can combine nicely

    1. Check if frame is rotated
    2. Use kmo_rotate to rotate the frame
    3. Update the headers
    4. Rename old file to _rotoffangle###
    5. Save file

    Note:
        Best to run this on sky sub files (_noS2.fits)

    Args:
        exposure (object): exposure object
        updateheader (bool): Update info in the header for the rotated frame
        keepsize (bool): extrapolate to keep the IFU the same size (for combining a mixture of cubes)


    """

    # Check if rotated
    if exposure.hdr['HIERARCH ESO OCS ROT OFFANGLE'] != 0:

        rotangle = exposure.hdr['HIERARCH ESO OCS ROT OFFANGLE']

        # Rotate
        if keepsize:
            extrapolate = '--extrapolate'
        else:
            extrapolate = ''

        kmo_rotate = 'esorex kmo_rotate --rotations="%f" %s %s' % (rotangle, extrapolate, exposure.filename)
        os.system(kmo_rotate)

        status, rotfile = subprocess.getstatusoutput("find . -maxdepth 1 -iname ROTATE.fits")

        # Rename files
        old_filename   = exposure.filename.replace('.fits', '_rotoffangle%.0f.fits' % (rotangle))
        rename_oldfile = 'mv %s %s' % (exposure.filename, old_filename)
        rename_newfile = 'mv %s %s' % (rotfile, exposure.filename)
        os.system(rename_oldfile)
        os.system(rename_newfile)

        if updateheader:

            newhdu = fits.open(exposure.filename)
            newhdr = newhdu[0].header
            newhdr['HIERARCH ESO OCS ROT OFFANGLE'] = 0.
            newhdr['HIERARCH ESO OCS ROT OFFSET']   = 0.

            if len(newhdu) > 3:
                for ifu in range(1, 25):
                    ext = newhdu['IFU.'+str(ifu)+'.DATA']
                    if len(ext.shape) > 0:
                        ext.header['CD1_1'] = ext.header['CDELT1']
                        ext.header['CD2_1'] = 0.
                        ext.header['CD1_2'] = 0.
                        ext.header['CD2_2'] = ext.header['CDELT2']
            else:
                for ext in range(1, 3):
                    newhdu[ext].header['CD1_1'] = newhdu[ext].header['CDELT1']
                    newhdu[ext].header['CD2_1'] = 0.
                    newhdu[ext].header['CD1_2'] = 0.
                    newhdu[ext].header['CD2_2'] = newhdu[ext].header['CDELT2']

            newhdu.writeto(exposure.filename, clobber=True)

        print(exposure.filename, 'rotation fixed')

    return


def flux_fix(exposure, clobber=True):
    """Fix Reflex exposure time bug

    Fix the bug in Reflex by dividing the calibrated flux from the standard stars
    by the DIT for science exposures

    Note:
        This is now obsolete as fixed in kmos pipeline v. > 1.4

    Args:
        exposure (object)
        clobber (bool): Overwrite exposure file


    """
    # Don't rewrite if we've already corrected this file
    if os.path.exists(exposure.filename_fluxfix) and clobber is False:
        print(exposure.filename_fluxfix, 'already exists, moving on...')

    else:
        count = 0
        for ifu in range(1, 25):
            ext = exposure.hdulist['IFU.'+str(ifu)+'.DATA']
            if ext.data is not None:
                # Count IFUs with data
                count = count + 1

                dat3D_fluxfix = ext.data / ext.header['EXPTIME']

                ext.data      = dat3D_fluxfix
                now           = datetime.datetime.now()
                ext.header['FLUX FIXED'] = str(now)

        print(count, 'IFUs with science data')

        exposure.hdulist.writeto(exposure.filename_fluxfix, clobber=clobber)
        print('Saved fixed flux to', exposure.filename_fluxfix)

    return
