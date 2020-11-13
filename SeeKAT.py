#!/usr/bin/env python
#Tiaan Bezuidenhout, 2020. For inquiries: bezmc93@gmail.com
#NB: REQUIRES Python 2

import argparse
import numpy as np
import math
import json
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
	--b 	Define the boresight to be different from the location of the
		brightest beam.
	'''

	parser.add_argument('-f', dest='file', 
				nargs = 1, 
				type = str, 
				help="Detections file",
				required=True)
	parser.add_argument('-c', dest='config', 
				nargs = 1, 
				type = str, 
				help="Configuration (json) file",
				required=False)
	parser.add_argument('-p',dest='psf',
				nargs=1,
				type=str,
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


def check_arguments(options):
	'''
	Checks that arguments make sense, and if so returns data read from input file. 
	'''

	data = np.genfromtxt(options.file[0],
			delimiter=' ',
			dtype=None,
			names=["RA","Dec","SN"],
			encoding="ascii")	
	best_cand = np.argsort(data["SN"])[-1]
	#print("Brightest beam: #"+str(best_cand))
	#print str(data["RA"][best_cand])
	#print str(data["Dec"][best_cand])
		
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
			return data,boresightCoord
		else:
			return data,boresightCoord
	else:
		ra, dec, options.overlap = read_from_json(options.config[0])
		b = SkyCoord(ra, dec, frame='icrs',unit=(u.hourangle, u.deg))
		boresightCoord = [b.ra.deg,b.dec.deg]

		if options.npairs[0] > n_comb:
			options.npairs[0] = n_comb
			return data,boresightCoord

		return data,boresightCoord


def read_from_json(config):
	'''
	Read tiling configuration from JSON file
	'''

	with open(config) as json_file:
		data = json.load(json_file)	
	ra = data['beams']['ca_target_request']['tilings'][0]['target'].split(',')[2]
	dec = data['beams']['ca_target_request']['tilings'][0]['target'].split(',')[3]
	overlap = str(data['beams']['ca_target_request']['tilings'][0]['overlap'])

	return ra, dec, overlap


def read_psf(psf):
	'''
	Converts PSF in fits format to a numpy array
	'''
	hdul = fits.open(psf)
	psf_ar = hdul[0].data

	#plt.imshow(psf_ar,origin='lower',cmap='gist_rainbow')
	#cbar = plt.colorbar()
	#cbar.ax.set_ylabel('Relative sensitivity')
	#plt.xlabel('Right Ascension (arcseconds)')
	#plt.ylabel('Declination (arcseconds)')
	#plt.show()


	#CLIPPING
	psf_ar[psf_ar<0.08] = 0

	return psf_ar


def convert_coords(c,psf,boresight,res):
	'''
	Converts coordinates from degrees to pixels
	'''
	### Build WCS 
	w = wcs.WCS(naxis=2)
	w.wcs.crpix = [0,0]
	step = res/3600.
	w.wcs.cdelt = np.array([-1*step, step])
	w.wcs.crval = boresight
	w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

	coordsDeg = []
	for i in range(0,len(c)):
		coordsDeg.append([c.ra.deg[i],c.dec.deg[i]])

	### Convert deg -> pix
	px = w.all_world2pix(coordsDeg,1)
	c.ra.px = px[:,0]
	c.dec.px = px[:,1]	

	boresightPx = w.all_world2pix([boresight],1)[0]
	### Calculate array size needed to fit all beams
	dist_right = np.absolute(max(c.ra.px) - boresightPx[0])
	dist_left = np.absolute(min(c.ra.px) - boresightPx[0])
	dist_up = np.absolute(max(c.dec.px) - boresightPx[1])
	dist_down = np.absolute(min(c.dec.px) - boresightPx[1])

	array_width = int(max(np.absolute(c.ra.px)*2) + psf.shape[0]+1)
	array_height = int(max(np.absolute(c.dec.px)*2) + psf.shape[0]+1)

	if array_width % 2 != 0: #making sure array has even dimensions
		array_width+=1
	if array_height % 2 != 0:
		array_height+=1

	#print("Width: "+str(array_width) + " px")
	#print("Height: "+str(array_height) + " px")

	### SHIFT coordinates so that boresight is at center of array
	c.ra.px += 0.5*array_width
	c.dec.px += 0.5*array_height 
	w.wcs.crpix = [0.5*array_width,0.5*array_height]

	#plt.scatter(c.ra.px,c.dec.px,color='lime')
	#plt.show()

	return c,w,array_width,array_height


def plot_known(w,known_coords):
	f, ax = plt.subplots()
	known_coords = known_coords.split(',')

	try:
		known_coords = [float(known_coords[0]),float(known_coords[1])]
	except:
		known_coords = SkyCoord(known_coords[0],known_coords[1], frame='icrs', unit=(u.hourangle, u.deg))
		known_coords = [known_coords.ra.deg,known_coords.dec.deg]


	known_px = w.all_world2pix([known_coords],1)
	
	plt.scatter(known_px[0,0],known_px[0,1],c='red',marker='x',s=200,zorder=999)	


def pix2deg(pix,w):
	coords_pix = []
	for i in range(0,len(pix[0])):
		coords_pix.append([pix[0][i],pix[1][i]])

	degs = w.all_pix2world(coords_pix,0)

	return degs

def get_best_pairs(data,npairs):
	sums = []
	for i in range (0,len(data)):
		for j in range (0,len(data)):
			if i !=j:
				sums.append(data["SN"][i]+data["SN"][j])

	threshold = sorted(sums)[-2*npairs]

	return threshold	


def make_plot(array_height,array_width,c,psf_ar,options,data):

	sum_threshold = get_best_pairs(data,options.npairs[0])

	full_ar = np.zeros((array_height,array_width))

	likelihood = np.zeros((array_height,array_width))

	for i in range(0,len(c)):
		#print i
		beam_ar = np.zeros((array_height,array_width))
		
		dec_start = int(c.dec.px[i])-int(psf_ar.shape[1])/2
		dec_end = int(c.dec.px[i])+int(psf_ar.shape[1])/2
		ra_start = int(c.ra.px[i])-int(psf_ar.shape[0])/2
		ra_end = int(c.ra.px[i])+int(psf_ar.shape[0])/2


		beam_ar[dec_start : dec_end,ra_start : ra_end] = psf_ar
		plt.contour(beam_ar,levels=[options.overlap],colors='white') # shows beam sizes

		full_ar = np.maximum(full_ar,beam_ar)


		for j in range(0,len(c)):

			stdout.write("\rComputing localisation curves for beam %d vs %d..." % (i,j))
			stdout.flush()
			if i!=j and data["SN"][i]+data["SN"][j] >= sum_threshold:
				
				plt.scatter(c.ra.px[i],c.dec.px[i],color='magenta')
				plt.scatter(c.ra.px[j],c.dec.px[j],color='magenta')
				comparison_ar = np.zeros((array_height,array_width))
		
				dec_start = int(c.dec.px[j])-int(psf_ar.shape[1])/2
				dec_end = int(c.dec.px[j])+int(psf_ar.shape[1])/2
				ra_start = int(c.ra.px[j])-int(psf_ar.shape[0])/2
				ra_end = int(c.ra.px[j])+int(psf_ar.shape[0])/2

				comparison_ar[dec_start : dec_end,
						ra_start : ra_end] = psf_ar

				plt.contour(comparison_ar,levels=[options.overlap],colors='white')

				beam_snr = data["SN"][i]
				comparison_snr = data["SN"][j]
				
				likelihood = localise(beam_snr,comparison_snr,beam_ar,comparison_ar,likelihood)


	plt.xlabel('RA ($^\circ$)')
	plt.ylabel('Dec ($^\circ$)')

	#plt.imshow(full_ar,origin='lower',cmap='inferno')
	#plt.show()	
	likelihood /= np.amax(likelihood)

	max_loc = np.where(likelihood==np.amax(likelihood))
	#max_loc =  np.transpose(max_loc)

	plt.imshow(likelihood,origin='lower',cmap='inferno')
	cbar = plt.colorbar()
	cbar.ax.set_ylabel('Localisation probability')
	plt.scatter(np.where(likelihood==np.amax(likelihood))[1],np.where(likelihood==np.amax(likelihood))[0],marker='v',c='cyan')
	plt.contour(likelihood,levels=[0.9],zorder=800,colors='cyan')
	plt.contour(likelihood,levels=[0.68],zorder=800,colors='lime')

	return max_loc

def get_ticks(array_width,array_height,w):
	num_ticks = min(array_width,array_height)/40
	ticks = [np.arange(0,min(array_width,array_height),min(array_width,array_height)/num_ticks), np.arange(0,min(array_width,array_height),min(array_width,array_height)/num_ticks)]

	labels = pix2deg(ticks,w)
	ra_deg= np.around(labels[:,0],4)
	dec_deg = np.around(labels[:,1],4)


	plt.xticks(ticks[0], ra_deg,rotation=40,fontsize=8)
	plt.yticks(ticks[1], dec_deg,fontsize=8)


def localise(beam_snr,comparison_snr,beam_ar,comparison_ar,likelihood):
	'''
	Plots contours where the ratio of the S/N detected in each 
	beam to the highest-S/N detection matches the ratio of 
	those beams' PSFs. 1-sigma errors are also drawn.
	'''
	ratio_ar = np.divide(beam_ar,comparison_ar)

	ratio_snr = beam_snr/comparison_snr		

	error = (1/beam_snr) + (1/comparison_snr)
				
	#plt.contour(ratio_ar, colors = 'white',
	#levels=[beam_snr/comparison_snr],zorder=800,origin=None)

	lower_bound = ratio_snr - error
	upper_bound = ratio_snr + error

	#ERROR CONTOURS
	#plt.contour(np.divide(beam_ar,comparison_ar),
	#linestyles='dashed',linewidths=0.7, colors = 'blue',
	#levels=[lower_bound], zorder=9000)

	#plt.contour(np.divide(beam_ar,comparison_ar),
	#linestyles='dashed',linewidths=0.7, colors = 'blue',
	#levels=[upper_bound], zorder=9000)

	gaussian = np.exp(-np.power(ratio_ar - ratio_snr, 2.) / (2 * np.power(error, 2.)))
	gaussian = np.nan_to_num(gaussian)	
	#plt.imshow(gaussian,cmap='inferno')
	#plt.show()

	likelihood += gaussian

	#plt.imshow(likelihood,cmap='inferno')
	#cbar = plt.colorbar()
	#cbar.ax.set_ylabel('Localisation probability')
	plt.xlabel('Right Ascension ($^\circ$)')
	plt.ylabel('Declination ($^\circ$)')
	return likelihood

def main():
	parser = argparse.ArgumentParser()
	options = parseOptions(parser)
	
	data,boresight = check_arguments(options)

	psf_ar = read_psf(options.psf[0])
	#psf_ar = psf_ar[0][0]	#FOR WSCLEAN PSF ONLY

	c = SkyCoord(data["RA"],data["Dec"], frame="icrs",unit=(u.hourangle,u.deg))

	c,w,array_width,array_height = convert_coords(c,psf_ar,boresight,options.res[0])
	
	if options.source:
		plot_known(w,options.source[0])

	ticks = get_ticks(array_width,array_height,w)

	max_loc = make_plot(array_height,array_width,c,psf_ar,options,data)
	max_deg = []
	

	print '\n'
	if len(max_loc) > 2:
		for m in max_loc:
			max_deg.append(pix2deg(m,w)[0])
		print '\nMaximum likelihood at coordinates:'
		print max_deg[0]


	plt.show()
main()

