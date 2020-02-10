#!/usr/bin/env python
#Tiaan Bezuidenhout, 2020. For inquiries: bezmc93@gmail.com
#NB: REQUIRES Python 2

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
from scipy.interpolate import griddata
from astropy.coordinates import SkyCoord
import astropy.units  as u
import argparse
import configparser
import datetime

import katpoint

from mosaic import PsfSim
from mosaic import beamforming
from mosaic import coordinate as coord


def parseOptions(parser):

	parser.add_argument('-f', dest='file', nargs = 1, 
		type = str, help="Candidate file",required=True)
	parser.add_argument('-c',dest='cfg',nargs = 1,
		type=str,help="Config file",required=True)
	parser.add_argument('-p',dest='psf',nargs=1,type=str,help="PSF file (optional)")
	options= parser.parse_args()
	return options



def parseConfig(cfgname):
	config = configparser.ConfigParser()
	config.read(cfgname)
	freq = float(config['DEFAULT']['FREQ'])
	coord = config['DEFAULT']['COORD'].split(',')
	time = config['DEFAULT']['TIME'].split(',')
	overlap = float(config['DEFAULT']['OVERLAP'])
	antennas = config['DEFAULT']['ANT'].split(',')

	if overlap > 1.0 or overlap < 0:
		print("The OVERLAP parameter must be between 0 and 1")
		exit()
	if not coord[1]:
		print("The COORD parameter must be in the format RA, Dec")
	if not time[1]:
		print("The TIME parameter must be in the format YYYY.MM.DD, hh:mm:ss")
	if not antennas[1]:
		print("The ANT parameter must be in the format of a comma-separated list")

	return freq,coord,time,overlap,antennas



def makeKatPointAntenna(antennaString):
	'''
	Takes in antenna number and returns katpoint Antenna object.
	'''

	antennaKat = []

	for antenna in antennaString:
		antkat = katpoint.Antenna(antenna)
		antennaKat.append(antkat)

	return antennaKat



def get_beamshape(freq,coord,time,overlap,antennas): 
	'''
	Given array parameters, writes coherent beam PSF to file, 
	and returns its shape.
	'''

	with open('antenna.csv', 'r') as antFile:
        	antennaCoords = antFile.readlines()

	antennaKat = makeKatPointAntenna([antennaCoords[int(ant)] for ant in antennas])

	psf = PsfSim(antennaKat, freq)

	boresight = katpoint.Target('boresight, radec, {}, {}'.format(coord[0], coord[1]))
	observeTime=datetime.datetime.strptime(time[0] + " " + time[1], '%Y.%m.%d %H:%M:%S.%f')
	beamShape = psf.get_beam_shape(boresight, observeTime)
	beamShape.psf.write_fits('psf.fits')
       
	el_width, el_height = beamShape.width_at_overlap(overlap)
		
	el_angle=-beamShape.angle
	
	return el_width,el_height,el_angle




def interpolate(psf_file,pulse,psf_px):
	'''Takes in a psf with a 20 x 20 pixel resolution, and interpolates
	to create a psf with resolution psf_px x psf_px.	
	'''

	hdul = fits.open(psf_file)
	data = hdul[0].data
	#plt.imshow(data)
	#plt.show()
	
	data*=max(pulse["SN"])
	psf_dim = data.shape
	points = np.zeros((psf_dim[0]*psf_dim[1],2))
	values = []
	n=0
	
	for i in range (0,psf_dim[0]):
		for j in range(0,psf_dim[1]):
			points[n] = [i, j]
			values.append(data[i,j])
			n+=1

	grid_x,grid_y = np.mgrid[0:19.0:19.0/psf_px,0:19.0:19.0/psf_px]
	model_prime = griddata(points,values,(grid_x,grid_y), method='cubic')
	
	#clipping model edges
	model_prime[model_prime<0.001*max(pulse["SN"])] = 0
	model_prime = np.rot90(model_prime,3)

	return model_prime



def pad_model(model_prime,p):
	'''
	Pads the PSF array for the brightest detection with zeroes
 	so that the other beams can be placed around it. 
	'''	

	n,m = model_prime.shape

	pad_min=int(np.floor(p/2))
	pad_max=int(np.ceil(float(p)/2))
	model_padded = np. zeros((p*n,p*m))
	model_padded[pad_min*n:pad_max*n,pad_min*m:pad_max*m]=model_prime

	return model_padded, pad_min, pad_max



def calc_padding(co,psf_px,psf_width,psf_height):
	'''
	This calculates how much padding to add to the primary PSF 
	based on how far apart the beams are.	 
	'''

	range_h = max(co.ra.deg) - min(co.ra.deg)
	range_v = max(co.dec.deg) - min(co.dec.deg)
	
	padding_fraction = (max(range_h,range_v) + 0.5*max(psf_width,psf_height))/max(psf_width,psf_height)
	padding_fraction = int(2*np.ceil(padding_fraction/2.0))
	if padding_fraction == 1:
		padding_fraction = 3
	if padding_fraction % 2 == 0:
		padding_fraction += 1
	print("Padding fraction = "+str(padding_fraction))

	return padding_fraction



def get_psf_size(model_prime,overlap,pulse):
	'''
	Determines the extent of the PSF (in pixels) at the level at 
	which the beams overlap. This is required to convert between
	pixels and real units. 
	'''

	model_cropped = model_prime.copy()
	model_cropped[model_prime<overlap*max(pulse["SN"])] = 0
	psf_height=max(np.nonzero(model_cropped)[0])- min(np.nonzero(model_cropped)[0])
	psf_width=max(np.nonzero(model_cropped)[1]) - min(np.nonzero(model_cropped)[1])
	
	return psf_height,psf_width



def convert_coords(pulse,psf_height,psf_width,el_height,el_width,psf_px):
	'''
	Converts coordinates from real units to pixels.
	'''

	#J1819_ra= 274.89233
	#J1819_dec = -14.96667
	#J1819_ra *= (psf_width/(2*el_width))
	#J1819_dec *= (psf_height/(2*el_height))

	c = SkyCoord(pulse["RA"],pulse["Dec"], frame="icrs",unit=(u.hourangle,u.deg))

	best_cand = np.argsort(pulse["SN"])[-1]
	best_coords=(c.ra.deg[best_cand],c.dec.deg[best_cand])

	c.ra.deg *= (psf_width/(2*el_width))  # coordinates now in pixels
	c.dec.deg *= (psf_height/(2*el_height))

	#J1819_ra -= c.ra.deg[best_cand]- (psf_px/2.0)*3 +psf_px/40.0	#When model_prime is rotated by -90 deg#
	#J1819_dec -= c.dec.deg[best_cand]- (psf_px/2.0)*3 -psf_px/40.0

	#plt.scatter(J1819_ra,J1819_dec,zorder=5545,c='yellow',marker='*')

	return c,best_cand,best_coords



def shift_coords(c,best_cand,psf_px,psf_height,psf_width,p):
	''''
	Shifts beams around in pixel space relative so that the brightest
	detection is in the centre, and the other beams are arranged around
	it.
	'''

	c.ra.deg -= c.ra.deg[best_cand]- (psf_px/2.0)*p +psf_px/40.0
	c.dec.deg -= c.dec.deg[best_cand]- (psf_px/2.0)*p -psf_px/40.0

	return c



def plot_coords(c,pulse):
	'''
	Plots detection coordinates with S/N colour scale
	'''

	sc=plt.scatter(c.ra.deg,c.dec.deg,c=pulse["SN"],zorder=47)
	#sc=plt.scatter(c.ra.deg,c.dec.deg,c='red',zorder=47)
	cb = plt.colorbar(sc)
	cb.set_label('S/N')
	ax = plt.gca()

	# Debugging tool: plot beam ellipses w/ el_width & el_height at correct positions.
	#for i in range(0,len(c.ra.deg)):
	#	e=ptch.Ellipse(xy=(c.ra.deg[i],c.dec.deg[i]),width=psf_width,
	#	height=psf_height,angle=el_angle,fill=False,color='white')
		
		#ax.add_artist(e)



def plot_psf(c,pulse,i,model_prime,model_padded,model_shifted,pad_min,pad_max,
	best_cand,n,m,psf_px,p,overlap):
	'''
	Shifts beam psf to appropriate position on model_shifted array, and 
	plots contour at overlap level.
	'''

	h_low = pad_min*n+int(c.ra.deg[i]-c.ra.deg[best_cand])
	h_high = pad_max*n+int(c.ra.deg[i]-c.ra.deg[best_cand])
	v_low = pad_min*m+int(c.dec.deg[i]-c.dec.deg[best_cand])
	v_high = pad_max*m+int(c.dec.deg[i]-c.dec.deg[best_cand])
	
	if h_low < 0 or h_high > p*psf_px or v_low < 0 or v_high > p*psf_px:
		print("Need more PADDING...")
		exit()

	model_shifted[v_low:v_high,h_low:h_high] += model_prime
		
	model = np.zeros(model_padded.shape)	

	model[v_low:v_high,h_low:h_high] = model_prime

	a = plt.contour(model,colors='black',levels=[overlap*max(pulse["SN"])])

	plt.imshow(model)
	
	psf = plt.imshow(model_shifted,vmin=8)

	return model_shifted,model



def localise(i,pulse,best_cand,model,model_padded):
	'''
	Plots contours where the ratio of the S/N detected in each 
	beam to the highest-S/N detection matches the ratio of 
	those beams' PSFs. 1-sigma errors are also drawn.
	'''

	if pulse["SN"][i]>0 and i!=best_cand:
		if len(pulse)>6:

			if pulse["SN"][i]>=pulse["SN"][np.argsort(pulse["SN"])[-6]]:
				error = (1/max(pulse["SN"]) + 1/pulse["SN"][i])
				
				psf = plt.contour(np.divide(model_padded,model),colors='red',
				levels=[max(pulse["SN"])/pulse["SN"][i]],zorder=800)

				#ERROR CONTOURS
				plt.contour(np.divide(model_padded,model),
				linestyles='dashed',linewidths=0.7,colors='magenta',
				levels=[max(pulse["SN"])/pulse["SN"][i]] - error, zorder=9000)

				plt.contour(np.divide(model_padded,model),
				linestyles='dashed',linewidths=0.7,colors='magenta',
				levels=[max(pulse["SN"])/pulse["SN"][i]] + error, zorder=9000)
		

		else:
			error = (1/max(pulse["SN"]) + 1/pulse["SN"][i])

			psf = plt.contour(np.divide(model_padded,model),colors='red',
			levels=[max(pulse["SN"])/pulse["SN"][i]], zorder=9000)
			
			#ERROR CONTOURS
			plt.contour(np.divide(model_padded,model),linestyles='dashed',linewidths=0.7,
			colors='magenta',levels=[max(pulse["SN"])/pulse["SN"][i]] - error, zorder=9000)

			plt.contour(np.divide(model_padded,model),linestyles='dashed',linewidths=0.7,
			colors='magenta',levels=[max(pulse["SN"])/pulse["SN"][i]] + error, zorder=9000)



def plot_beams(c,pulse,model_prime,model_padded,pad_min,pad_max,best_cand,best_coords,
	overlap,el_height,el_width,el_angle,psf_height,psf_width,p,filename,psf_px):
	'''
	Main plotting function.
	'''
		
	n,m = model_prime.shape
	model_shifted=np.zeros((p*n,p*m))

	for i in range (0,len(pulse)):
		print(i)

		model_shifted,model = plot_psf(c,pulse,i,model_prime,model_padded,
		model_shifted,pad_min,pad_max,best_cand,n,m,psf_px,p,overlap)
		
		localise(i,pulse,best_cand,model,model_padded)	

	locs = np.arange(0,psf_px*p,psf_px*p/5)
#	labels = [round((float(item)-psf_px/2.0*p -psf_px/40.0)*(2*el_width/psf_width)+best_coords[0],4) for item in locs]
	labels = [round((float(item)-psf_px/2.0*p +psf_px/40.0)*(2*el_width/psf_width)+
		best_coords[0],4) for item in locs] # When model_prime is rotated -90 deg 

	plt.xticks(locs, labels,rotation=60)
	labels = [round((float(item)-psf_px/2.0*p -psf_px/40.0)*(2*el_height/psf_height)+
		best_coords[1],5) for item in locs] 

	plt.yticks(locs, labels)
	plt.grid(linestyle='-')
	plt.autoscale()
	plt.gca().invert_yaxis()
	plt.title(filename)

	plt.xlabel("Right Ascension (deg)") 
	plt.ylabel("Declination (deg)")
	plt.tight_layout()
	plt.savefig("Contours_"+filename+".png",dpi=300)
	plt.show()
	plt.close()



def main():
	psf_px = 1000
	parser = argparse.ArgumentParser()
	options = parseOptions(parser)
	
	if options.psf != None:
		psf_file = options.psf[0]
	else:
		psf_file = 'psf.fits'
	
	pulse = np.genfromtxt(options.file[0],delimiter=' ',dtype=None,
		names=["RA","Dec","SN"],encoding="ascii")
		

	freq,coord,time,overlap,antennas = parseConfig(options.cfg[0])

	el_width,el_height,el_angle = get_beamshape(freq,coord,time,overlap,antennas)
	
	model_prime = interpolate(psf_file, pulse,psf_px)
	
	psf_height,psf_width = get_psf_size(model_prime,overlap,pulse)

	print("Converting coordinates...")
	c,best_cand,best_coords = convert_coords(pulse,psf_height,psf_width,el_height,el_width,psf_px)

	p = calc_padding(c,psf_px,psf_width,psf_height)

	c = shift_coords(c,best_cand,psf_px,psf_height,psf_width,p)

	model_padded,pad_min,pad_max = pad_model(model_prime,p)
	
	print("Plotting...")
	plot_coords(c,pulse)
	
	plot_beams(c,pulse,model_prime,model_padded,pad_min,
		pad_max,best_cand,best_coords,overlap,
		el_height,el_width,el_angle,psf_height,
		psf_width,p,options.file[0],psf_px)


main()






