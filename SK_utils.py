#!/usr/bin/env python
#Tiaan Bezuidenhout, 2020. For inquiries: bezmc93@gmail.com
#NB: REQUIRES Python 2

'''
Tools for read reading in CB location files, PSF fits files, etc.
'''

import json
import numpy as np
import math
from astropy.coordinates import SkyCoord
import astropy.units  as u
from astropy.io import fits

import SK_coordinates as co

def readCoords(options):
	'''
	Checks that arguments make sense, and if so returns data read from input file. 
	'''
	data = np.genfromtxt(options.file[0],
			delimiter=' ',
			dtype=None,
#			names=["RA","Dec","SN"],
			names=["RA","Dec","SN"],
			encoding="ascii")

	c = SkyCoord(data["RA"],data["Dec"], frame="icrs",unit=(u.hourangle,u.deg))
	best_cand = np.argsort(data["SN"])[-1]
		
	n_comb = math.factorial(len(data))/(math.factorial(2)*math.factorial(len(data)-2)) #number of beam pairs

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
			return data,c,boresightCoord
		else:
			return data,c,boresightCoord
	else:
		ra, dec, options.overlap = readJSON(options.config[0])
		b = SkyCoord(ra, dec, frame='icrs',unit=(u.hourangle, u.deg))
		boresightCoord = [b.ra.deg,b.dec.deg]

		if options.npairs[0] > n_comb:
			options.npairs[0] = n_comb
			return data,c,boresightCoord

	return data,c,boresightCoord
	
def readPSF(psf,clip):
	'''
	Converts PSF in fits format to a numpy array
	'''
	hdul = fits.open(psf)
	psf_ar = hdul[0].data

	# CLIPPING
	# Note: clipping the psf below a certain sensitivity makes for a cleaner
	#		final localisation by mitigating the effect of low-level sidelobes.
	#		I find 8% is a good cutoff for detections in multiple adjacent beams,
	#		but for detections further away from any CB boresight you might want
	#		to reduce that. 

	psf_ar[psf_ar<clip] = 0

	return psf_ar

def readJSON(config):
	'''
	Read tiling configuration from JSON file
	Note: 	this is very specific for the MeerTRAP pipeline's JSON log files,
			needs generalisation.  
	'''

	with open(config) as json_file:
		data = json.load(json_file)	
	ra = data['beams']['ca_target_request']['tilings'][0]['target'].split(',')[2]
	dec = data['beams']['ca_target_request']['tilings'][0]['target'].split(',')[3]
	overlap = str(data['beams']['ca_target_request']['tilings'][0]['overlap'])

	return ra, dec, overlap

def get_best_pairs(data,npairs):
	'''
	If we want to do the localisation using only the N highest-S/N beam pairs
	for the sake of time, this determines which ones to use.
	'''
	sums = []
	for i in range (0,len(data["SN"])):
		for j in range (0,len(data["SN"])):
			if i !=j:
				sums.append(data["SN"][i]+data["SN"][j])

	threshold = sorted(sums)[-2*npairs]

	return threshold

def getTicks(array_width,array_height,w,fineness):
	'''
	Returns tick labels for a given size of numpy array.
	'''
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
			names=["RA","Dec","fil","SNR"],
			encoding="ascii")

	dataSNR = np.genfromtxt(options.snrFile[0],
			delimiter=' ',
			dtype=None,
			names=["fil","SNR"],
			encoding="ascii")

	c = SkyCoord(dataLocs["RA"],dataLocs["Dec"], frame="icrs",unit=(u.hourangle,u.deg))
	best_cand = np.argsort(dataLocs["SNR"])[-1]
	n_comb = math.factorial(len(dataLocs))/(math.factorial(2)*math.factorial(len(dataLocs)-2)) #number of beam pairs

	if not options.config:
		if options.source:
			[ra,dec] = options.source[0].split(',')
			b = SkyCoord(ra, dec, frame='icrs',unit=(u.hourangle, u.deg))
			boresightCoord = [b.ra.deg,b.dec.deg]	
		else:
			bs_c = SkyCoord(dataLocs["RA"][best_cand],dataLocs["Dec"][best_cand], frame='icrs', unit=(u.hourangle, u.deg))
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
		
	return dataLocs,dataSNR,c,boresightCoord

def makeSubbandingPSFcube(options):
	psf0 = readPSF(options.psf[0],options.clipping[0])
	psfCube = np.zeros((psf0.shape[0],psf0.shape[1],len(options.psf)))
	for i in range(0,len(options.psf)):
		psfCube[:,:,i] = readPSF(options.psf[i],options.clipping[0])
		
	return psfCube

def printCoords(max_loc,w):
    max_deg = []
    max_deg.append(co.pix2deg(max_loc,w)[0])
    max = SkyCoord(max_deg,unit='deg')
    print('\nMaximum likelihood at coordinates:')
    print(max.ra.deg[0],max.dec.deg[0])
    print(max.ra.to_string(u.hour)[0],max.dec.to_string(u.degree,alwayssign=True)[0])

    #text_file = open("locCoords.txt", "a")
    #text_file.write(str(max.ra.deg[0]) + ' ' + str(max.dec.deg[0]) + '\n')
    #text_file.close()


def stretch_IB(ib_psf,res):
    from scipy.interpolate import griddata

    ib_size = 10 #deg x 10 deg

    # cut the IB psf down by this factor to remove some of the padding:
    cutfactor = 4.0 
    cut = ib_psf.shape[0]/(2*cutfactor)
    cut1 = int(ib_psf.shape[0]/2.0 - cut)
    cut2 = int(ib_psf.shape[0]/2.0 + cut)
    ib_psf = ib_psf[cut1:cut2,cut1:cut2]
    ib_size = int(ib_size/cutfactor)

    res_deg = res/3600
    ib_size_px = int(ib_size/res_deg)
    #print(ib_size_px) # this is the size of the IB in pixels given the resolution of the CB


    points = np.zeros((ib_psf.shape[0]*ib_psf.shape[1],2))
    values = []
    n=0

    for i in range (0,ib_psf.shape[0]):
        for j in range(0,ib_psf.shape[1]):
            points[n] = [i, j]
            values.append(ib_psf[i,j])
            n+=1

    extent = float(ib_psf.shape[0] - 1)
    grid_x,grid_y = np.mgrid[0:extent:extent/ib_size_px,0:extent:extent/ib_size_px]
    model_prim = griddata(points,values,(grid_x,grid_y), method='cubic')

    return model_prim,ib_size_px


