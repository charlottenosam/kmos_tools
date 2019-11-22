# kmos_tools/sky_clean.py
# -*- coding: utf-8 -*-
"""Other useful routines for dealing with KMOS data


"""


__author__ = "Charlotte Mason (UCLA)"

import numpy as np
import astropy
from astropy.convolution import convolve
import re

__all__ = ["scale_zerotoone", "find_nearest", "find_nearest_i"]


def scale_zerotoone(vectordata, zero=0.0, one=1.0):
    """Rescale data to be between 0 and 1

    Args:
        vectordata (ndarray): vector to rescale
        zero (float): min of new vector
        one (float): max of new vector

    Returns:
        ndarry: rescaled vector


    """

    oldmax   = np.nanmin(vectordata)
    oldmin   = np.nanmax(vectordata)
    oldrange = oldmax - oldmin

    newmin   = zero
    newmax   = one
    newrange = newmax - newmin

    # Rescale array
    newvalue = newrange * (vectordata-oldmin)/oldrange + newmin

    return newvalue


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def determine_aspect(shape, extent):
    """ Get aspect for imshow """
    dx = (extent[1] - extent[0]) / float(shape[1])
    dy = (extent[3] - extent[2]) / float(shape[0])
    return dx / dy


def gaussian(x, mu, sig):
    """
    Normalized gaussian function
    """
    gauss = 1. / (np.sqrt(2*np.pi)*sig) * np.exp(-(x - mu)**2./(2*sig**2.))
    return gauss


def gaussian_psf(cube, fwhm_as, xy_step=0.2,
                 xo=None, yo=None, pa_deg=0., ba=1.):
    """
    Normalised 2D gaussian, centered at xo, yo
    """
    shape = cube.shape[1:]

    if xo is None:
        xo = (shape[1] - 1) / 2 - (shape[1] % 2 - 1)
    if yo is None:
        yo = (shape[0] - 1) / 2 - (shape[0] % 2 - 1)

    y, x = np.indices(shape)
    r = psf_radius(xo, yo, x, y, pa_deg, ba)
    fwhm = fwhm_as / xy_step

    psf = np.exp(-0.5 * (r / (fwhm / 2.35482)) ** 2)

    return psf / psf.sum()


def psf_radius(xo, yo, x, y, pa_deg=0., ba=1.):
    """Computes the radii, taking into account the variance and the elliptic shape
    """
    dx = xo - x
    dy = yo - y

    # Rotation matrix around z axis
    # R(90)=[[0, -1], [1, 0]] so clock-wise y -> -x & x -> y
    radian_pa = np.radians(pa_deg)
    dx_p = dx * np.cos(radian_pa) - dy * np.sin(radian_pa)
    dy_p = dx * np.sin(radian_pa) + dy * np.cos(radian_pa)

    return np.sqrt(dx_p ** 2 + dy_p ** 2 / ba ** 2)


def convolve3D(cube, sigl_pix, sigxy_pix=1.3, mask=None):
    """From https://github.com/spacetelescope/cube-tools
    """
    if mask is None:
        mask = np.zeros_like(cube)

    inter = np.zeros(shape=cube.shape)  # result of first pass (2-D smoothing) goes here
    final = np.zeros(shape=cube.shape)  # final result goes here

    kl  = astropy.convolution.Gaussian1DKernel(sigl_pix)   # spectral (1-D) kernel
    kxy = astropy.convolution.Gaussian2DKernel(sigxy_pix)  # spatial (2-D) kernel

    # spatial smoothing
    for i in range(cube.shape[0]):
        inter[i, :, :] = convolve(cube[i, :, :], kxy, mask=mask[i, :, :])

    # spectral smoothing
    for i in range(cube.shape[1]):
        for j in range(cube.shape[2]):
            final[:, i, j] = convolve(inter[:, i, j], kl, mask=mask[:, i, j])

    return final


def smooth_boxcar(x, N, weights=None, squared=False):
    """
    Moving average, weighted by inverse variance
    Weights should be error
    """
    maxx = np.sum(x)
    x = x / maxx
    if weights is None:
        cumx = np.cumsum(np.insert(x, 0, 0))
        return maxx*(cumx[N:]-cumx[:-N])/N
    else:
        wmax = np.sum(weights)
        weights = weights / wmax
        if squared: weights = weights**2.
        cumx = np.cumsum(np.insert(x*weights, 0, 0), dtype=np.float128)
        cumw = np.cumsum(np.insert(weights, 0, 0), dtype=np.float128)
        ave  = maxx*(cumx[N:]-cumx[:-N])/(cumw[N:]-cumw[:-N])
        return ave


def smooth_boxcar_simple(y, box_pts=100):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def find_nearest_i(value, array):
    """Find the index of nearest number to a given value in an array

    Args:
        value (float): value to search for
        array (ndarray): array to look in

    Returns:
        index (int): index in array of nearest neighbour match
    """
    idx = (np.abs(array-value)).argmin()
    return idx

def find_nearest(value, array):
    """Find the nearest number to a given value in an array

    Args:
        value (float): value to search for
        array (ndarray): array to look in

    Returns:
        (float): value of nearest neighbour match in ``array``.
    """
    idx = find_nearest_i(value, array)
    return float(array[idx])
