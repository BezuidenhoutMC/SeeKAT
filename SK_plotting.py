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
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

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

def likelihoodPlot(ax,loglikelihood):
    '''
    Creates the localisation plot
    '''
    likelihood = np.exp(loglikelihood - np.nanmax(loglikelihood))
    likelihood /= likelihood.sum()

    plt.imshow(likelihood,origin='lower',cmap='inferno')
    
    scalebar = AnchoredSizeBar(ax.transData,
                           10, '10 arcseconds', 'upper left', 
                           pad=1.0,
                           color='black',
                           frameon=True,
                           size_vertical=0.2)

    ax.add_artist(scalebar)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Localisation likelihood')
    #cbar.ax.set_ylabel('Relative intensity')


    plt.scatter(np.where(likelihood==np.amax(likelihood))[1],np.where(likelihood==np.amax(likelihood))[0],marker='v',c='cyan')

    #plt.contour(likelihood,levels=[1-np.std(likelihood)],zorder=800,colors='cyan')
    #plt.contour(likelihood,levels=[1-2*np.std(likelihood)],zorder=800,colors='lime')

    ## Calculating the interval values in 2D
    likelihood_flat_sorted = np.sort(likelihood, axis=None)
    likelihood_flat_sorted_cumsum = np.cumsum(likelihood_flat_sorted)
    
    if len(np.nonzero(likelihood_flat_sorted_cumsum > (1-0.6827))[0]) != 0:
        ind_1sigma = np.nonzero(likelihood_flat_sorted_cumsum > (1-0.6827))[0][0]
        val_1sigma = likelihood_flat_sorted[ind_1sigma]
        plt.contour(likelihood,levels=[val_1sigma],zorder=800,colors='cyan')

    if len(np.nonzero(likelihood_flat_sorted_cumsum > (1-0.9545))[0]) != 0:
        ind_2sigma = np.nonzero(likelihood_flat_sorted_cumsum > (1-0.9545))[0][0]
        val_2sigma = likelihood_flat_sorted[ind_2sigma]
        plt.contour(likelihood,levels=[val_2sigma],zorder=800,colors='lime')

    #if len(np.nonzero(likelihood_flat_sorted_cumsum > (1-0.9973))[0]) != 0:
        #ind_3sigma = np.nonzero(likelihood_flat_sorted_cumsum > (1-0.9973))[0][0]
        #val_3sigma = likelihood_flat_sorted[ind_3sigma]
        #plt.contour(likelihood,levels=[val_3sigma],zorder=800,colors='green')

    #max_loc = np.where(likelihood==np.amax(likelihood))
    ## When displaying a 2D array, the last index is the last axis, thus we need to flip things here.
    #max_loc = [max_loc[1],max_loc[0]]


    plt.xlabel('RA ($^\circ$)')
    plt.ylabel('Dec ($^\circ$)')