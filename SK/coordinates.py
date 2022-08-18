#!/usr/bin/env python
# Tiaan Bezuidenhout, 2020. For inquiries: bezmc93@gmail.com
# NB: REQUIRES Python 3

"""
Tools for converting between real and pixel coordinates.
"""

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy.wcs as wcs


def buildWCS(boresight, step, resolution):
    """
    Builds a World Coordinate System around the observation boresight
    """

    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [0, 0]
    step = resolution / 3600.0
    w.wcs.cdelt = np.array([-1 * step, step])
    w.wcs.crval = boresight
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    return w


def calcArraySize(c, boresightPx, psf):
    """
    Calculates the size of numpy array needed to fit all the CBs simultaneously
    c must be a SkyCoord object
    """
    dist_right = np.absolute(max(c.ra.px) - boresightPx[0])
    dist_left = np.absolute(min(c.ra.px) - boresightPx[0])
    dist_up = np.absolute(max(c.dec.px) - boresightPx[1])
    dist_down = np.absolute(min(c.dec.px) - boresightPx[1])

    array_width = int(max(np.absolute(c.ra.px) * 2) + psf.shape[0] + 1)
    array_height = int(max(np.absolute(c.dec.px) * 2) + psf.shape[0] + 1)
    
    if array_width % 2 != 0: # making sure array has even dimensions
        array_width += 1
    if array_height % 2 != 0:
        array_height += 1

    return array_width, array_height


def deg2pix(c, psf, boresight, res):
    """
    Converts coordinates from degrees to pixels
    c must be a SkyCoord object and psf a numpy array.
    """

    step = res / 3600.
    w = buildWCS(boresight, step, res)
    coordsDeg = []
    for i in range(0, len(c)):
        coordsDeg.append([c.ra.deg[i], c.dec.deg[i]])

    ### Convert deg -> pix
    px = w.all_world2pix(coordsDeg, 1)
    c.ra.px = px[:, 0]
    c.dec.px = px[:, 1]
    ### Calculate array size needed to fit all beams
    boresightPx = w.all_world2pix([boresight], 1)[0]

    array_width,array_height = calcArraySize(c, boresightPx, psf)

    ### SHIFT coordinates so that boresight is at center of array
    c.ra.px += 0.5 * array_width
    c.dec.px += 0.5 * array_height 
    w.wcs.crpix = [0.5 * array_width, 0.5 * array_height]

    return c, w,array_width, array_height

def pix2deg(pix, w):
    """
    Converts pixel numbers to coordinates according to a pre-defined WCS.
    """
    
    coords_pix = []
    for i in range(0, len(pix[0])):
        coords_pix.append([pix[0][i], pix[1][i]])
    degs = w.all_pix2world(coords_pix, 0)
    
    return degs

def convert_coords_IB(c, psf, boresight, res):
    '''
    Converts coordinates from degrees to pixels
    '''
    
    ### Build WCS 
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [0, 0]
    step = res / 3600.0
    w.wcs.cdelt = np.array([-1 * step, step])
    w.wcs.crval = boresight
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    coordsDeg = []
    for i in range(0, len(c)):
        coordsDeg.append([c.ra.deg[i], c.dec.deg[i]])

    ### Convert deg -> pix
    px = w.all_world2pix(coordsDeg, 1)
    c.ra.px = px[:, 0]
    c.dec.px = px[:, 1]

    boresightPx = w.all_world2pix([boresight],1)[0]

    ### Calculate array size needed to fit all beams
    array_width = psf.shape[0]
    array_height = psf.shape[1]

    ### SHIFT coordinates so that boresight is at center of array
    c.ra.px += 0.5 * array_width
    c.dec.px += 0.5 * array_height 
    w.wcs.crpix = [0.5 * array_width, 0.5 * array_height]

    return c,w