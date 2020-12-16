#!/usr/bin/env python
#Tiaan Bezuidenhout, 2020. For inquiries: bezmc93@gmail.com
#NB: REQUIRES Python 2

import SK_utils as ut
import SK_coordinates as co
import SK_plotting as Splot
import SeeKAT_v5 as SK

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
	-f	Input file with each line a different CB detection.
		Should have 3 columns: RA (h:m:s), Dec (d:m:s), S/N
	-p	PSF of a CB in fits format
	--o	Fractional sensitivity level at which CBs are tiled to overlap
	--r	Resolution of PSF in units of arcseconds per pixel
	--n	Number of beams to consider when creating overlap contours. Will
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
	parser.add_argument('--s', dest='source',
				nargs = 1,
				type=str,
				help="Draws given coordinate location (degrees) on localisation plot",
				required = False)


	options= parser.parse_args()

	return options



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    options = parseOptions(parser)
    
    #data,c,boresight = ut.readCoords(options)
    
    dataLocs,dataSNR,c,boresight = ut.readSubbandingFiles(options)
    
    psfCube = ut.makeSubbandingPSFcube(options.psf)
    
    c,w,array_width,array_height = co.deg2pix(c,psfCube[:,:,0],boresight,options.res[0])
    
    f, ax = plt.subplots()
    
    if options.source:
        Splot.plot_known(w,options.source[0])

    Splot.make_ticks(array_width,array_height,w,fineness=100)
    
    fullbandLikelihood = np.zeros((array_height,array_width))
    for subband in range(0,psfCube.shape[2]):
        print "\n---Localising in subband %d/%d---" % (subband+1,psfCube.shape[2])
        subbandSNRs = []
        for line in dataSNR:
            if line["fil"][-2:] == "."+str(subband):
                subbandSNRs.append(line["SNR"])
        subbandData = {"SN":np.asarray(subbandSNRs)}

        likelihood = SK.make_plot(array_height,array_width,c,psfCube[:,:,subband],options,subbandData)
       #plt.imshow(likelihood)
        #plt.show()
        
        fullbandLikelihood+=likelihood*subbandData["SN"][subband]
        fullbandLikelihood /= np.amax(fullbandLikelihood)

    
    Splot.likelihoodPlot(fullbandLikelihood)
    max_deg = []
    max_loc = np.where(likelihood==np.amax(likelihood))
    #max_loc =  np.transpose(max_loc)
    print '\n'
    if len(max_loc) > 2:
        for m in max_loc:
            max_deg.append(pix2deg(m,w)[0])
        print '\nMaximum likelihood at coordinates:'
        print max_deg[0]

    plt.show()
