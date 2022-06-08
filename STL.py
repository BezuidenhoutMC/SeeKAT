#!/usr/bin/env python
from SK import utils as ut
from SK import coordinates as co
from SK import plotting as Splot
import SeeKAT as SK

import argparse
import numpy as np
import math
import matplotlib.pyplot as plt
from sys import stdout
from astropy.io import fits

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

def make_plot(array_height,array_width,c,psf_ar,options,data,header):
    full_ar = np.empty((array_height,array_width))

    likelihood = np.empty((array_height,array_width))
    likelihood += 1.

    det_idx = 0
    det_ar = np.empty((array_height,array_width))

    dec_start = int(np.round(c.dec.px[det_idx]))-int(psf_ar.shape[1]/2)
    dec_end = int(np.round(c.dec.px[det_idx]))+int(psf_ar.shape[1]/2)
    ra_start = int(np.round(c.ra.px[det_idx]))-int(psf_ar.shape[0]/2)
    ra_end = int(np.round(c.ra.px[det_idx]))+int(psf_ar.shape[0]/2)

    det_ar[dec_start : dec_end,ra_start : ra_end] = psf_ar
    full_ar = det_ar
    c1 = plt.contour(det_ar,levels=[options.overlap],colors='black',linewidths=1.5) # shows beam sizes

    for con in c1.collections[0].get_paths():
        xs = []
        ys = []
        v = con.vertices
        #print '-------'
        #print v

        xs = v[:,0]
        ys = v[:,1]

        v_degs = co.pix2deg((xs,ys),w)
        #print v_degs
        line = 'polygon'
        for q in range(0,len(v_degs)):
            line += ' ' + str(v_degs[q,0])
            line += ' ' + str(v_degs[q,1])
        header.append(line)
        #print c.ra.to_string(u.hour)[i].replace('h',':').replace('m',':').replace('s',''),c.dec.to_string(u.degree,alwayssign=True)[i].replace('d',':').replace('m',':').replace('s',''), beam_ar[565,554]*10.0/0.814798480193

    for i in range(0,len(c)):
        #print i

        stdout.write("\rComputing localisation curves for beam %d/%d..." % (i+1,len(c)))
        stdout.flush()
        if i!=det_idx:
                plt.scatter(c.ra.px[i],c.dec.px[i],color='black',s=2.)

                comparison_ar = np.zeros((array_height,array_width))
                dec_start = int(np.round(c.dec.px[i]))-int(psf_ar.shape[1]/2)
                dec_end = int(np.round(c.dec.px[i]))+int(psf_ar.shape[1]/2)
                ra_start = int(np.round(c.ra.px[i]))-int(psf_ar.shape[0]/2)
                ra_end = int(np.round(c.ra.px[i]))+int(psf_ar.shape[0]/2)

                comparison_ar[dec_start : dec_end,
                    ra_start : ra_end] = psf_ar
                full_ar	= np.maximum(full_ar,comparison_ar)

                c2 = plt.contour(comparison_ar,levels=[options.overlap],colors='black',linewidths=0.5,linestyles='dashed')

                det_snr = data["SN"][det_idx]
                comparison_snr = data["SN"][i]

                likelihood = localise(det_snr,comparison_snr,det_ar,comparison_ar,likelihood)
                
    likelihood[det_ar == 0] = 0.0

    likelihood /= np.amax(likelihood)

    likelihood*=full_ar


    return likelihood,header


def localise(beam_snr,comparison_snr,beam_ar,comparison_ar,likelihood):
    '''
    Plots contours where the ratio of the S/N detected in each 
    beam to the highest-S/N detection matches the ratio of 
    those beams' PSFs. 1-sigma errors are also drawn.
    '''
    det_thresh = 8.0
    ratio_ar = np.divide(beam_ar,comparison_ar)

    likelihood[ratio_ar <= beam_snr/det_thresh] *= 0.
    return likelihood

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    options = parseOptions(parser)
    
    data,c,boresight = ut.readCoords(options)
    psf_ar = ut.readPSF(options.psf[0],clip=0.0)
    

    c,w,array_width,array_height = co.deg2pix(c,psf_ar,boresight,options.res[0])
    
    f, ax = plt.subplots(figsize=(10,10))
    
    if options.source:
        Splot.plot_known(w,options.source[0])

    Splot.make_ticks(ax,array_width,array_height,w,fineness=50) #fineness=x places ticks every x pixels

    header = ['# Region file format: DS9 version 4.1',
          ('global color=blue dashlist=8 3 '
           'width=2 font="helvetica 10 '
           'normal roman" select=1 highlite=1 '
           'dash=0 fixed=0 edit=1 move=1 '
           'delete=1 include=1 source=1'),
          'fk5']


    loglikelihood,header = make_plot(array_height,array_width,c,psf_ar,options,data,header)
    #likelihood[0:int(0.45*array_height),0:int(array_width)] = 0.
    #likelihood[array_height-int(0.45*array_height):array_height,0:int(array_width)] = 0.
    #likelihood[0:int(array_height),0:int(0.45*array_width)] = 0.
    #likelihood[0:int(array_height),array_width-int(0.45*array_width):int(array_width)] = 0.

    #likelihood[array_width-int(0.1*array_width):array_width,array_height-int(0.1*array_height):array_height] = 0.
    #likelihood[0:200,0:1330] = 0
    #likelihood[:] = 0

    # likelihood /= likelihood.sum()

    Splot.likelihoodPlot(f,ax,loglikelihood,options)

    print('Writing to fits...\n')
    h = w.to_header()
    hdu = fits.PrimaryHDU(header=h)
    hdu.data = likelihood
    hdu.update_header
    hdu.writeto('localisation.fits',overwrite=True, output_verify='fix')

    plt.savefig(options.file[0]+'_localisation.png',dpi=300)
    plt.show()
