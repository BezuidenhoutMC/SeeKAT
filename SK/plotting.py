#!/usr/bin/env python
# Tiaan Bezuidenhout, 2020. For inquiries: bezmc93@gmail.com
# NB: REQUIRES Python 3

"""
Plotting tools.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

import SK.utils as ut
import SK.coordinates as co


def plot_known(w, known_coords):
    """
    Makes a cross on the localisation plot where a known source is located.
    w must be a WCS object pre-set with SK_coordinates.buildWCS
    """

    known_coords = known_coords.split(",")

    try:
        known_coords = [float(known_coords[0]), float(known_coords[1])]
    except:
        known_coords = SkyCoord(
            known_coords[0], known_coords[1], frame="icrs", unit=(u.hourangle, u.deg)
        )
        known_coords = [known_coords.ra.deg, known_coords.dec.deg]

    known_px = w.all_world2pix([known_coords], 1)

    plt.scatter(known_px[0, 0], known_px[0, 1], c="cyan", marker="x", s=20, zorder=999)


def make_ticks(ax, array_width, array_height, w, fineness):
    """
    Adds ticks and labels in sky coordinates to the plot.
    """

    ticks = ut.getTicks(array_width, array_height, w, fineness)
    labels = co.pix2deg(ticks, w)
    ra_deg = np.around(labels[:, 0], 4)
    dec_deg = np.around(labels[:, 1], 4)
    ax.set_xticks(ticks[0])  # ,rotation=40,fontsize=8)
    ax.set_xticklabels(ra_deg, rotation=40)
    ax.set_yticks(ticks[1])  # , dec_deg)#,fontsize=8)
    ax.set_yticklabels(dec_deg)


def likelihoodPlot(f, ax, loglikelihood, options):
    """
    Creates the localisation plot
    """
    # np.save('CB_localisation.npy',loglikelihood)

    likelihood = np.exp(loglikelihood - np.nanmax(loglikelihood))
    likelihood = np.nan_to_num(likelihood)

    likelihood /= likelihood.sum()

    likelihood[likelihood < 0.000000001] = 0.0
    plt.imshow(likelihood, origin="lower", cmap="binary", aspect="auto")
    plt.subplots_adjust(left=0.1, right=0.85, bottom=0.1, top=0.85)

    ## Making adjustable scale bar
    if options.sb[0] != 0:
        scalebar = AnchoredSizeBar(
            ax.transData,
            int(options.sb[0] / options.res[0]),
            "%d arcsec" % (options.sb[0]),
            "upper left",
            pad=1.0,
            color="black",
            frameon=True,
            size_vertical=0.2,
        )

        ax.add_artist(scalebar)

    # cbar = plt.colorbar()

    plt.scatter(
        np.where(likelihood == np.amax(likelihood))[1],
        np.where(likelihood == np.amax(likelihood))[0],
        marker="x",
        c="red",
    )

    ## Calculating the interval values in 2D
    likelihood_flat_sorted = np.sort(likelihood, axis=None)
    likelihood_flat_sorted_cumsum = np.cumsum(likelihood_flat_sorted)

    if len(np.nonzero(likelihood_flat_sorted_cumsum > (1 - 0.6827))[0]) != 0:
        ind_1sigma = np.nonzero(likelihood_flat_sorted_cumsum > (1 - 0.6827))[0][0]
        val_1sigma = likelihood_flat_sorted[ind_1sigma]
        CS1 = plt.contour(likelihood, levels=[val_1sigma], zorder=800, colors="red")
    else:
        val_1sigma = np.nan

    if len(np.nonzero(likelihood_flat_sorted_cumsum > (1 - 0.9545))[0]) != 0:
        ind_2sigma = np.nonzero(likelihood_flat_sorted_cumsum > (1 - 0.9545))[0][0]
        val_2sigma = likelihood_flat_sorted[ind_2sigma]
        CS2 = plt.contour(
            likelihood, levels=[val_2sigma], zorder=800, linestyles="--", colors="red"
        )
    else:
        val_2sigma = np.nan

    if len(np.nonzero(likelihood_flat_sorted_cumsum > (1 - 0.9973))[0]) != 0:
        ind_3sigma = np.nonzero(likelihood_flat_sorted_cumsum > (1 - 0.9973))[0][0]
        val_3sigma = likelihood_flat_sorted[ind_3sigma]
        # CS3 = plt.contour(likelihood,levels=[val_3sigma],zorder=800,colors='lightgray')
    else:
        val_3sigma = np.nan

    lines = [CS1.collections[0], CS2.collections[0]]
    labels = ["1-$\sigma$", "2-$\sigma$"]
    plt.legend(lines, labels, loc="upper right", borderpad=1.0, edgecolor="black")

    ## Making axis histograms
    ax_histx = f.add_axes([0.1, 0.855, 0.75, 0.1], sharex=ax)
    ax_histx.plot(np.sum(likelihood, axis=0), color="black")
    plt.setp(ax_histx.get_xticklabels(), visible=False)

    # ax_histx.fill_between(x=range(0,likelihood.shape[1]),y1=-0.1*max(np.sum(likelihood,axis=0)),y2=np.sum(likelihood,axis=0),where=np.amax(likelihood,axis=0)>=val_3sigma,color='magenta')
    # ax_histx.fill_between(x=range(0,likelihood.shape[1]),y1=-0.1*max(np.sum(likelihood,axis=0)),y2=np.sum(likelihood,axis=0),where=np.amax(likelihood,axis=0)>=val_2sigma,color='lime')
    # ax_histx.fill_between(x=range(0,likelihood.shape[1]),y1=-0.1*max(np.sum(likelihood,axis=0)),y2=np.sum(likelihood,axis=0),where=np.amax(likelihood,axis=0)>=val_1sigma, color='cyan')

    ax_histy = f.add_axes([0.855, 0.1, 0.1, 0.75], sharey=ax)
    ax_histy.plot(
        np.sum(likelihood, axis=1), range(0, likelihood.shape[0]), color="black"
    )
    ax_histy.set_title("Likelihood", size=10)
    plt.setp(ax_histy.get_yticklabels(), visible=False)

    # ax_histy.fill_betweenx(y=range(0,likelihood.shape[0]),x1=-0.1*max(np.sum(likelihood,axis=1)),x2=np.sum(likelihood,axis=1),where=np.amax(likelihood,axis=1)>=val_3sigma, color='magenta')
    # ax_histy.fill_betweenx(y=range(0,likelihood.shape[0]),x1=-0.1*max(np.sum(likelihood,axis=1)),x2=np.sum(likelihood,axis=1),where=np.amax(likelihood,axis=1)>=val_2sigma, color='lime')
    # ax_histy.fill_betweenx(y=range(0,likelihood.shape[0]),x1=-0.1*max(np.sum(likelihood,axis=1)),x2=np.sum(likelihood,axis=1),where=np.amax(likelihood,axis=1)>=val_1sigma, color='cyan')

    ax.set_xlabel("RA ($^\circ$)")
    ax.set_ylabel("Dec ($^\circ$)")
