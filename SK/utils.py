#!/usr/bin/env python
# Tiaan Bezuidenhout, 2020. For inquiries: bezmc93@gmail.com
# NB: REQUIRES Python 3

"""
Tools for read reading in CB location files, PSF fits files, etc.
"""

import json
import numpy as np
import math
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from scipy.interpolate import griddata

import SK.coordinates as co


def readCoords(options):
    """
    Checks that arguments make sense, and if so returns data read from input file. 
    """

    data = np.genfromtxt(options.file[0],
            delimiter=' ',
            dtype=None,
            names=["RA", "Dec", "SN"],
            encoding="ascii")

    # Sorting beams by S/N, so that Gaussian assumption holds better. 
    # You generally want LOW/HIGH S/N beam pairs.
    data = data[np.argsort(data["SN"])[::-1]] 

    c = SkyCoord(data["RA"],data["Dec"], frame="icrs",unit=(u.hourangle, u.deg))
    best_cand = np.argsort(data["SN"])[-1]

    # Calculate number of beam pairs
    n_comb = math.factorial(len(data))/(math.factorial(2)*
    						math.factorial(len(data)-2)) 

    if not options.config:
        if options.source:
            [ra,dec] = options.source[0].split(',')
            b = SkyCoord(ra, dec, frame='icrs',unit=(u.hourangle, u.deg))
            boresightCoord = [b.ra.deg,b.dec.deg]
        else:
            bs_c = SkyCoord(data["RA"][best_cand],data["Dec"][best_cand], frame='icrs', unit=(u.hourangle, u.deg))
            boresightCoord = [bs_c.ra.deg,bs_c.dec.deg]

        if options.overlap > 1.0 or options.overlap < 0:
            print("The OVERLAP parameter must be between 0 and 1")
            exit()

        elif options.npairs[0] > n_comb:
            options.npairs[0] = n_comb
            return data, c, boresightCoord
        else:
            return data, c, boresightCoord
    else:
        ra, dec, options.overlap = readJSON(options.config[0])
        b = SkyCoord(ra, dec, frame='icrs',unit=(u.hourangle, u.deg))
        boresightCoord = [b.ra.deg, b.dec.deg]

        if options.npairs[0] > n_comb:
            options.npairs[0] = n_comb
            return data, c, boresightCoord
        else:
            return data, c, boresightCoord


def readPSF(psf,clip):
    """
    Converts PSF in fits format to a numpy array
    """

    hdul = fits.open(psf)
    psf_ar = hdul[0].data

    # CLIPPING
    # Note: clipping the psf below a certain sensitivity makes for a cleaner
    #       final localisation by mitigating the effect of low-level sidelobes.
    #       I find 8% is a good cutoff for detections in multiple adjacent beams,
    #       but for detections further away from any CB boresight you might want
    #       to reduce that. 

    psf_ar[psf_ar < clip] = 0.0

    return psf_ar


def norm_likelihood(loglikelihood):
    likelihood = np.exp(loglikelihood - np.nanmax(loglikelihood))
    likelihood = np.nan_to_num(likelihood)
    likelihood /= likelihood.sum()
    likelihood[likelihood<0.000000001] = 0.0

    return likelihood


def CDF(s):
    # Cumulative distribution function for 2D Gaussian
    return  1 - np.exp(-0.5 * s ** 2)


def calc_error(likelihood, sigma):
    likelihood_flat_sorted = np.sort(likelihood, axis=None)
    likelihood_flat_sorted_index = np.argsort(likelihood, axis=None)
    likelihood_flat_sorted_cumsum = np.cumsum(likelihood.flatten()[likelihood_flat_sorted_index])

    if len(np.nonzero(likelihood_flat_sorted_cumsum > (1 - CDF(sigma)))[0]) != 0:
        # Index where cum sum goes above sigma percentage
        index_sigma_above = np.nonzero(likelihood_flat_sorted_cumsum > 
                        (1 - CDF(sigma)))[0]

        # Minimum likelihood included in error
        likelihood_index_sigma_sorted = likelihood_flat_sorted_index[index_sigma_above]

        likelihood_sigma_level = likelihood_flat_sorted[index_sigma_above[0]]

        likelihood_index_sigma_original = np.unravel_index(likelihood_index_sigma_sorted, 
                                                            likelihood.shape)
        y_lower_index, y_higher_index = np.sort(likelihood_index_sigma_original[0])[[np.s_[0], 
                                                np.s_[-1]]]
        x_lower_index, x_higher_index = np.sort(likelihood_index_sigma_original[1])[[np.s_[0], 
                                                np.s_[-1]]]
        max_likehood_index = np.unravel_index(np.argmax(likelihood), likelihood.shape)

        x_lower = x_lower_index, max_likehood_index[0]
        x_higher = x_higher_index, max_likehood_index[0]
        y_lower = max_likehood_index[1], y_lower_index
        y_higher = max_likehood_index[1], y_higher_index

        return likelihood_sigma_level, [x_lower, x_higher, y_lower, y_higher]

    else:
        print('Could not find error')
        exit()


def to_HMS_second(coord):
        return coord.hms.h*60*60 + coord.hms.m*60 + coord.hms.s

def to_DMS_second(coord):
        return coord.dms.d*60*60 + coord.dms.m*60*60 + coord.dms.s


def printCoords(max_loc, w):
    if len(max_loc) == 2:
        max_loc = (max_loc[1],max_loc[0])
        max_deg = []
        max_deg.append(co.pix2deg(max_loc, w)[0])
        max = SkyCoord(max_deg, unit='deg')
        print('\nMaximum likelihood at coordinates:')
        print('{:.3f}, {:.3f} deg'.format(max.ra.deg[0], max.dec.deg[0]))

        print(max.ra.to_string(u.hour)[0] + ', ' + 
                max.dec.to_string(u.degree, alwayssign=True)[0])
        print('')


def printError(max_loc, w, error, sigma):
    if len(max_loc) == 2:
        max_loc = (max_loc[1],max_loc[0])

        max_deg = []
        max_deg.append(co.pix2deg(max_loc, w)[0])
        center_eq = SkyCoord(max_deg, unit='deg')

        error_eq_deg =  w.all_pix2world(error, 0)
        error_eq = SkyCoord(error_eq_deg, unit='deg')

        ra_lower = center_eq.ra - error_eq[0].ra
        ra_higher = error_eq[1].ra - center_eq.ra
        dec_lower = center_eq.dec - error_eq[2].dec
        dec_higher = error_eq[3].dec - center_eq.dec
    
        print("| RA error (s): -{:.3f}, +{:.3f}".format(
                np.abs(to_HMS_second(ra_lower)[0]),
                np.abs(to_HMS_second(ra_higher)[0])))
        print("| DEC error (arcs): -{:.3f}, +{:.3f}".format(
                np.abs(to_DMS_second(dec_lower)[0]),
                np.abs(to_DMS_second(dec_higher)[0])))
        print('')


def write2fits(w, likelihood):
    header = w.to_header()
    hdu = fits.PrimaryHDU(header=header)
    hdu.data = likelihood
    hdu.update_header
    fitsname = 'Likelihood.fits'
    hdu.writeto(fitsname, overwrite=True)


def readJSON(config):
    """
    Read tiling configuration from JSON file
    Note: this is very specific for the MeerTRAP pipeline's JSON log files,
          needs generalisation.  
    """

    with open(config) as json_file:
        data = json.load(json_file)	
    ra = data['beams']['ca_target_request']['tilings'][0]['target'].split(',')[2]
    dec = data['beams']['ca_target_request']['tilings'][0]['target'].split(',')[3]
    overlap = str(data['beams']['ca_target_request']['tilings'][0]['overlap'])

    return ra, dec, overlap

def get_best_pairs(data,npairs):
    """
    If we want to do the localisation using only the N highest-S/N beam pairs
    for the sake of time, this determines which ones to use.
    """

    sums = []
    for i in range (0,len(data["SN"])):
        for j in range (0,len(data["SN"])):
            if i !=j:
                sums.append(data["SN"][i] + data["SN"][j])

    threshold = sorted(sums)[-2*npairs]

    return threshold


def getTicks(array_width,array_height,w,fineness):
    """
    Returns tick labels for a given size of numpy array.
    """

    try:
        num_ticks = min(array_width,array_height)/fineness
        ticks = [np.arange(0,min(array_width,array_height),min(array_width,array_height)/num_ticks), np.arange(0,min(array_width,array_height),min(array_width,array_height)/num_ticks)]

    except:
        num_ticks = 10
        ticks = [np.arange(0,min(array_width,array_height),min(array_width,array_height)/num_ticks), np.arange(0,min(array_width,array_height),min(array_width,array_height)/num_ticks)]

    return ticks


def readSubbandingFiles(options):
    dataLocs = np.genfromtxt(options.locFile[0],
            delimiter=' ',
            dtype=None,
            names=["RA", "Dec", "fil", "SNR"],
            encoding="ascii")

    dataSNR = np.genfromtxt(options.snrFile[0],
            delimiter=' ',
            dtype=None,
            names=["fil", "SNR"],
            encoding="ascii")

    c = SkyCoord(dataLocs["RA"],dataLocs["Dec"], frame="icrs",unit=(u.hourangle,u.deg))
    best_cand = np.argsort(dataLocs["SNR"])[-1]
    n_comb = math.factorial(len(dataLocs)) / (math.factorial(2)*
            math.factorial(len(dataLocs) - 2)) #number of beam pairs

    if not options.config:
        if options.source:
            [ra,dec] = options.source[0].split(',')
            b = SkyCoord(ra, dec, frame='icrs',unit=(u.hourangle, u.deg))
            boresightCoord = [b.ra.deg, b.dec.deg]	
        else:
            bs_c = SkyCoord(dataLocs["RA"][best_cand],
                            dataLocs["Dec"][best_cand], 
                            frame='icrs', unit=(u.hourangle, u.deg))
            boresightCoord = [bs_c.ra.deg,bs_c.dec.deg]

        if options.overlap > 1.0 or options.overlap < 0:
            print("The OVERLAP parameter must be between 0 and 1")
            exit()

        elif options.npairs[0] > n_comb:
            options.npairs[0] = n_comb
            return dataLocs,dataSNR,c,boresightCoord
        else:
            return dataLocs,dataSNR,c,boresightCoord
    else:
        ra, dec, options.overlap = readJSON(options.config[0])
        b = SkyCoord(ra, dec, frame='icrs',unit=(u.hourangle, u.deg))
        boresightCoord = [b.ra.deg,b.dec.deg]

        if options.npairs[0] > n_comb:
            options.npairs[0] = n_comb
            return dataLocs,dataSNR,c,boresightCoord
        else:
            return data, c, boresightCoord

    return dataLocs,dataSNR,c,boresightCoord


def makeSubbandingPSFcube(options):
    psf0 = readPSF(options.psf[0], options.clipping[0])
    psfCube = np.zeros((psf0.shape[0], psf0.shape[1], len(options.psf)))
    for i in np.arange(0, len(options.psf)):
        psfCube[:,:,i] = readPSF(options.psf[i], options.clipping[0])

    return psfCube


def autozoom(ax, c, options):
    ax.set_xlim(min(c.ra.px) - int(15/options.res[0]), 
                max(c.ra.px) + int(15/options.res[0]))
    ax.set_ylim(min(c.dec.px) - int(15/options.res[0]), 
                max(c.dec.px) + int(15/options.res[0]))


def stretch_IB(ib_psf, res):
    ib_size = 10 #deg x 10 deg

    # cut the IB psf down by this factor to remove some of the padding:
    cutfactor = 4.0 
    cut = ib_psf.shape[0] / (2 * cutfactor)
    cut1 = int(ib_psf.shape[0] / 2.0 - cut)
    cut2 = int(ib_psf.shape[0] / 2.0 + cut)
    ib_psf = ib_psf[cut1:cut2, cut1:cut2]
    ib_size = int(ib_size / cutfactor)

    res_deg = res / 3600.0
    ib_size_px = int(ib_size/res_deg)

    points = np.zeros((ib_psf.shape[0]*ib_psf.shape[1], 2))
    values = []
    n = 0

    for i in range (0, ib_psf.shape[0]):
        for j in np.arange(0, ib_psf.shape[1]):
            points[n] = [i, j]
            values.append(ib_psf[i, j])
            n += 1

    extent = float(ib_psf.shape[0] - 1)
    grid_x,grid_y = np.mgrid[0:extent:extent / ib_size_px, 0:extent:extent / ib_size_px]
    model_prim = griddata(points, values, (grid_x, grid_y), method='cubic')

    return model_prim, ib_size_px