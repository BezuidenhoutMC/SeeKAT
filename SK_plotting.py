#!/usr/bin/env python
#Tiaan Bezuidenhout, 2020. For inquiries: bezmc93@gmail.com
#NB: REQUIRES Python 2

'''
Plotting tools.
'''

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units  as u

import SK_utils as ut
import SK_coordinates as co


def plot_known(w,known_coords):
	'''
	Makes a cross on the localisation plot where a known source is located.
	w must be a WCS object pre-set with SK_coordinates.buildWCS
	'''

	known_coords = known_coords.split(',')

	try:
		known_coords = [float(known_coords[0]),float(known_coords[1])]
	except:
		known_coords = SkyCoord(known_coords[0],known_coords[1], frame='icrs', unit=(u.hourangle, u.deg))
		known_coords = [known_coords.ra.deg,known_coords.dec.deg]


	known_px = w.all_world2pix([known_coords],1)
	
	plt.scatter(known_px[0,0],known_px[0,1],c='red',marker='x',s=200,zorder=999)

def make_ticks(array_width,array_height,w,fineness):
	'''
	Adds ticks and labels in sky coordinates to the plot.
	'''

	ticks = ut.getTicks(array_width,array_height,w,fineness)
	labels = co.pix2deg(ticks,w)
	ra_deg= np.around(labels[:,0],4)
	dec_deg = np.around(labels[:,1],4)
	plt.xticks(ticks[0], ra_deg,rotation=40,fontsize=8)
	plt.yticks(ticks[1], dec_deg,fontsize=8)

def likelihoodPlot(ax,likelihood):
	'''
	Creates the localisation plot
	'''

	plt.imshow(likelihood,origin='lower',cmap='inferno')

	cbar = plt.colorbar()
	cbar.ax.set_ylabel('Localisation probability')
	plt.scatter(np.where(likelihood==np.amax(likelihood))[1],np.where(likelihood==np.amax(likelihood))[0],marker='v',c='cyan')
#	plt.contour(likelihood,levels=[0.9],zorder=800,colors='cyan')
#	plt.contour(likelihood,levels=[0.68],zorder=800,colors='lime')
	#plt.contour(likelihood,levels=[0.9],zorder=800,colors='cyan')
	plt.contour(likelihood,levels=[0.32],zorder=800,colors='lime')

	plt.xlabel('RA ($^\circ$)')
	plt.ylabel('Dec ($^\circ$)')
	
	#ax.set_xlim(450,670)
	#ax.set_ylim(450,670)
