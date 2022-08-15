#!/usr/bin/env python
#Tiaan Bezuidenhout, 2020. For inquiries: bezmc93@gmail.com
#NB: REQUIRES Python 3

import SK.utils as ut
import SK.coordinates as co
import SK.plotting as Splot
import SeeKAT as SK

import argparse
import numpy as np
import math
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units  as u
import astropy.wcs as wcs
from sys import stdout
import matplotlib.patches as patches


def parseOptions(parser):
    '''Options:
    -f    Input file
    --o    Fractional sensitivity level at which CBs are tiled to overlap
    --r    Resolution of PSF in units of arcseconds per pixel
    --n    Number of beams to consider when creating overlap contours. Will
        pick the specified number of beams with the highest S/N values.
    '''

    parser.add_argument('-l', dest='locFile', 
                nargs = 1, 
                type = str, 
                required=True)
    parser.add_argument('-c', dest='config', 
                nargs = 1, 
                type = str, 
                help="Configuration (json) file",
                required=False)
    parser.add_argument('-s', dest='snrFile', 
                nargs = 1, 
                type = str, 
                required=True)
    parser.add_argument('-p',dest='psf',
                                type=str,
                                nargs = '+', 
                                help="PSF file",
                                required=True)
    parser.add_argument('--o', dest='overlap',
                type = float,
                help = "Fractional sensitivity level at which the coherent beams overlap",
                default = 0.25,
                required = False)
    parser.add_argument('--r', dest='res',
                nargs = 1,
                type=float,
                help="Distance in arcseconds represented by one pixel of the PSF",
                default = 1,
                required = True)
    parser.add_argument('--n',dest='npairs',
                nargs = 1,
                type = int,
                default = [1000000])
    parser.add_argument('--nsig',dest='nsig',
                nargs = 1,
                type = int,
                help='Draws uncertainty contours up to this number of standard deviations.',
                default = [2])
    parser.add_argument('--s', dest='source',
                nargs = 1,
                type=str,
                help="Draws given coordinate location (degrees) on localisation plot",
                required = False)
    parser.add_argument('--scalebar', dest='sb',
                        nargs = 1,
                        type = float,
                        help = "Length of scale bar on the localisation plot in arcseconds. Set to 0 to omit altogether",
                        default = [10],
                        required = False)
    parser.add_argument('--ticks', dest='tickspacing',
                        nargs = 1,
                        type = float,
                        help = "Sets the number of pixels between ticks on the localisation plot",
                        default = [50],
                        required = False)
    parser.add_argument('--clip', dest='clipping',
                        nargs = 1,
                        type = float,
                        help = "Sets values of the PSF below this number to zero. Helps minimise the influence of low-level sidelobes",
                        default = [0.08],
                        required = False)    
    parser.add_argument('--fits', dest='fitsOut',
                        help = "Outputs .fits file of localisation region",
                        action = 'store_true')   

    options= parser.parse_args()

    return options



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    options = parseOptions(parser)
    
    #data,c,boresight = ut.readCoords(options)
    
    dataLocs,dataSNR,c,boresight = ut.readSubbandingFiles(options)
    
    psfCube = ut.makeSubbandingPSFcube(options)
    
    c,w,array_width,array_height = co.deg2pix(c,psfCube[:,:,0],boresight,options.res[0])
    
    f, ax = plt.subplots(figsize=(10,10))
    
    if options.source:
        Splot.plot_known(w,options.source[0])

    Splot.make_ticks(ax,array_width,array_height,w,fineness=options.tickspacing[0])
    
    fullbandLogLikelihood = np.zeros((array_height,array_width))
    for subband in range(0,psfCube.shape[2]):
        print("\n---Localising in subband %d/%d---" % (subband+1,psfCube.shape[2]))
        subbandSNRs = []
        for line in dataSNR:
            if line["fil"][-2:] == "."+str(subband):
                subbandSNRs.append(line["SNR"])
        subbandData = {"SN":np.asarray(subbandSNRs)}
        loglikelihood = SK.make_map(array_height,array_width,c,psfCube[:,:,subband],options,subbandData)
       #plt.imshow(likelihood)
        #plt.show()

        fullbandLogLikelihood+=loglikelihood#*np.max(subbandData["SN"])
        #fullbandLogLikelihood /= np.amax(fullbandLikelihood)
  
    Splot.likelihoodPlot(f,ax, w, fullbandLogLikelihood,options)
    max_deg = []
    max_loc = np.where(fullbandLogLikelihood==np.nanmax(fullbandLogLikelihood))
    if len(max_loc) == 2:
        max_loc = (max_loc[1],max_loc[0])
        ut.printCoords(max_loc,w)
    else:
        print('Multiple equally possible locations')

    plt.show()
