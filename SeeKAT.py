#!/usr/bin/env python
#Tiaan Bezuidenhout, 2020. For inquiries: bezmc93@gmail.com
#NB: REQUIRES Python 2

import argparse
import numpy as np
import matplotlib.pyplot as plt
from sys import stdout

import SK_utils as ut
import SK_coordinates as co
import SK_plotting as Splot


def parseOptions(parser):
	'''Options:
	-f	Input file with each line a different CB detection.
		Should have 3 columns: RA (h:m:s), Dec (d:m:s), S/N
	-p	PSF of a CB in fits format
	--o	Fractional sensitivity level at which CBs are tiled to overlap
	--r	Resolution of PSF in units of arcseconds per pixel
	--n	Number of beams to consider when creating overlap contours. Will
		pick the specified number of beams with the highest S/N values.
	--s Draws known coordinates onto the plot for comparison.
	--scalebar Sets the length of the scale bar on the plot in arcseconds.
	--ticks Sets the spacing of ticks on the localisation plot.
	--clip Sets level below which CB PSF is set equal to zero.
	--zoom Automatically zooms in on the TABs.
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
	parser.add_argument('--zoom', dest='autozoom',
						help = "Automatically zooms the localisation plot in on the TABs",
						action = 'store_true')	
	options= parser.parse_args()

	return options


def make_plot(array_height,array_width,c,psf_ar,options,data):
    
    sum_threshold = ut.get_best_pairs(data,int(options.npairs[0]))   # Only beam pairs with S/N summing to above this
                                                                # number will be used for localisation

    full_ar = np.zeros((array_height,array_width))

    loglikelihood = np.zeros((array_height,array_width))

    for i in range(0,len(c)):
        #print i
        beam_ar = np.zeros((array_height,array_width))

        dec_start = int(c.dec.px[i])-int(psf_ar.shape[1]/2)
        dec_end = int(c.dec.px[i])+int(psf_ar.shape[1]/2)
        ra_start = int(c.ra.px[i])-int(psf_ar.shape[0]/2)
        ra_end = int(c.ra.px[i])+int(psf_ar.shape[0]/2)

        beam_ar[dec_start : dec_end,ra_start : ra_end] = psf_ar
        plt.contour(beam_ar,levels=[options.overlap],colors='white',linewidths=0.5,linestyles='dashed') # shows beam sizes

        full_ar = np.maximum(full_ar,beam_ar)

        for j in range(0,len(c)):
            if i<j and data["SN"][i]+data["SN"][j] >= sum_threshold:
                stdout.write("\rComputing localisation curves for beam %d vs %d/%d..." % (j+1,i+1,len(c)))
                stdout.flush()
                plt.scatter(c.ra.px,c.dec.px,color='white',s=0.2)

                comparison_ar = np.zeros((array_height,array_width))

                dec_start = int(c.dec.px[j])-int(psf_ar.shape[1]/2)
                dec_end = int(c.dec.px[j])+int(psf_ar.shape[1]/2)
                ra_start = int(c.ra.px[j])-int(psf_ar.shape[0]/2)
                ra_end = int(c.ra.px[j])+int(psf_ar.shape[0]/2)

                comparison_ar[dec_start : dec_end, ra_start : ra_end] = psf_ar

                plt.contour(comparison_ar,levels=[options.overlap],colors='white',linewidths=0.5)
                plt.contour(beam_ar,levels=[options.overlap],colors='white',linewidths=0.5)

                beam_snr = data["SN"][i]
                comparison_snr = data["SN"][j]

                loglikelihood = localise(beam_snr,comparison_snr,beam_ar,comparison_ar,loglikelihood)


    #loglikelihood /= np.amax(loglikelihood)
    #plt.imshow(full_ar,origin='lower',cmap='inferno')

    return loglikelihood


def localise(beam_snr,comparison_snr,beam_ar,comparison_ar,loglikelihood):
    '''
    Plots contours where the ratio of the S/N detected in each 
    beam to the highest-S/N detection matches the ratio of 
    those beams' PSFs. 1-sigma errors are also drawn.
    '''

    ratio_ar = np.divide(beam_ar,comparison_ar)

    ratio_snr = beam_snr/comparison_snr		

    #error = (1/beam_snr) + (1/comparison_snr)   # Old error formula, Gaussian approximation
    error = (1.0 / comparison_snr) + (beam_snr / comparison_snr**2)  # Full error formula, https://en.wikipedia.org/wiki/Ratio_distribution

    lower_bound = ratio_snr - error
    upper_bound = ratio_snr + error
    
    gaussian = np.exp(-np.power(ratio_ar - ratio_snr, 2.) / (2 * np.power(error, 2.)))
    gaussian /= 2.0*np.pi*error
    gaussian = np.nan_to_num(gaussian)
    gaussian = np.log(gaussian)

    loglikelihood += gaussian

    return loglikelihood

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    options = parseOptions(parser)
    
    data,c,boresight = ut.readCoords(options)
    
    psf_ar = ut.readPSF(options.psf[0],options.clipping[0])
    
    c,w,array_width,array_height = co.deg2pix(c,psf_ar,boresight,options.res[0])
    
    f, ax = plt.subplots()
    
    if options.source:
        Splot.plot_known(w,options.source[0])

    Splot.make_ticks(array_width,array_height,w,fineness=options.tickspacing[0])
    
    loglikelihood = make_plot(array_height,array_width,c,psf_ar,options,data)
    #np.save("CB_localisation.npy",loglikelihood)

    Splot.likelihoodPlot(ax,loglikelihood,options)
    max_deg = []
    max_loc = np.where(loglikelihood==np.amax(loglikelihood))
    
    if len(max_loc) == 2:
        max_loc = (max_loc[1],max_loc[0])
        ut.printCoords(max_loc,w)
    else:
        print('Multiple equally possible locations')
    
    if options.autozoom == True:
        ax.set_xlim(min(c.ra.px) - int(15/options.res[0]) ,max(c.ra.px) + int(15/options.res[0]))
        ax.set_ylim(min(c.dec.px) -int(15/options.res[0]) ,max(c.dec.px) + int(15/options.res[0]))

    plt.savefig(options.file[0]+'.png',dpi=300)
    plt.show()
