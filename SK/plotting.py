#!/usr/bin/env python
# Tiaan Bezuidenhout, 2020. For inquiries: bezmc93@gmail.com
# NB: REQUIRES Python 3

"""
Plotting tools.
"""
import matplotlib.pyplot as plt

# Make text larger for readability on graphs
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['font.size'] = 22
# Set same fonts as my LaTeX document
plt.rcParams['font.family'] = 'STIXGeneral' 
plt.rcParams['mathtext.fontset'] = 'stix'
# Other plotting params
plt.rcParams['figure.figsize'] = [10,10]

import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

import SK.utils as ut
import SK.coordinates as co


def plot_known(w,known_coords):
    """
    Makes a cross on the localisation plot where a known source is located.
    w must be a WCS object pre-set with SK_coordinates.buildWCS
    """
    known_coords = known_coords.split(',')

    try:
        known_coords = [float(known_coords[0]), float(known_coords[1])]
    except:
        known_coords = SkyCoord(known_coords[0], known_coords[1], frame='icrs', 
                                unit=(u.hourangle, u.deg))
        known_coords = [known_coords.ra.deg, known_coords.dec.deg]

    known_px = w.all_world2pix([known_coords], 1)

    plt.scatter(known_px[0, 0],known_px[0, 1], c='cyan', marker='x', s=20, zorder=999)


def make_ticks(ax, array_width, array_height, w, fineness):
    """
    Adds ticks and labels in sky coordinates to the plot.
    """

    ticks = ut.getTicks(array_width, array_height, w, fineness)
    labels = co.pix2deg(ticks, w)
    ra_deg= np.around(labels[:, 0], 4)
    dec_deg = np.around(labels[:, 1], 4)
    ax.set_xticks(ticks[0])
    ax.set_xticklabels(ra_deg)
    ax.set_yticks(ticks[1])
    ax.set_yticklabels(dec_deg)


def likelihoodPlot(f, ax, w, loglikelihood, options):
    """
    Creates the localisation plot
    """
    
    likelihood = ut.norm_likelihood(loglikelihood)

    import matplotlib.colors
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#FFFFFF","#1E88E5"])

    plt.imshow(likelihood, origin='lower', cmap=cmap, aspect="auto")
    plt.subplots_adjust(left=0.17, right=0.85, bottom=0.1, top=0.85)

    plt.scatter(np.where(likelihood==np.amax(likelihood))[1],
                np.where(likelihood==np.amax(likelihood))[0], 
                marker='x', c='red')

    ## Making adjustable scale bar 
    fontprops = fm.FontProperties(size=18)   
    if options.sb[0] != 0:
        scalebar = AnchoredSizeBar(ax.transData,
                           int(options.sb[0]/options.res[0]), 
                           '%d arcsec' % (options.sb[0]), 
                           'upper left', 
                           pad=1.0,
                           color='#1d3557',
                           frameon=True,
                           size_vertical=0.2,
                           fontproperties=fontprops)

        ax.add_artist(scalebar)
    
    # Printing location of maximum likelihood
    max_loc = np.where(loglikelihood==np.nanmax(loglikelihood))
    ut.printCoords(max_loc, w)
    # print(max_loc)
    if len(max_loc[0]) == 2:
        plt.axhline(max_loc[0], lw=1.2, c='#a8dadc', ls='--')
        plt.axvline(max_loc[1], lw=1.2, c='#a8dadc', ls='--')

    ## Calculating the interval values in 2D
    print('Error levels:')
    sigs = np.asarray(np.arange(1,options.nsig[0]+1))

    for s in sigs:
        level, error = ut.calc_error(likelihood, s)

        if s == 1:
            ls = '-'
            lc = '#e63946'
            lw = 3

        else:
            ls = '--'
            lc = '#E65476'
            lw = 2

        plt.contour(likelihood, levels=[level], zorder=800,
                        colors=lc, linestyles=ls, linewidths=lw)

        print('---- %i sigma ----' % (s))
        ut.printError(max_loc, w, error, s)

    ## Making axis histograms
    ax_histx = f.add_axes([0.17, 0.855, 0.68, 0.1], sharex=ax)
    ax_histx.plot(np.sum(likelihood, axis=0),color='black')
    plt.setp(ax_histx.get_xticklabels(), visible=False)
    if len(max_loc[0]) == 2:
        ax_histx.vlines(max_loc[1], 0, np.max(np.sum(likelihood, axis=0)) ,
                     lw=1.2, colors='#a8dadc', ls='--')
    
    ax_histy = f.add_axes([0.855, 0.1, 0.1, 0.75], sharey=ax)
    ax_histy.plot(np.sum(likelihood, axis=1),np.arange(0,likelihood.shape[0]), 
                    color='black')
    ax_histy.set_title('Likelihood', size=15)
    plt.setp(ax_histy.get_yticklabels(), visible=False)
    if len(max_loc[0]) == 2:
        ax_histy.hlines(max_loc[0], 0, np.max(np.sum(likelihood, axis=1)) ,
                     lw=1.2, colors='#a8dadc', ls='--')

    ax.set_xlabel("RA (deg)")
    ax.set_ylabel("Dec (deg)")


    
