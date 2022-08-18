import argparse
import numpy as np
import math
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

import sys
#sys.path.append('/raid/tbez/SeeKAT/') #REMOVE

import SK.plotting as Splot
import SK.utils as ut
import SK.coordinates as co
import SeeKAT as SK

def parseOptions(parser):
    '''Options:
    -f    Input files with each line a different CB detection.
        Should have 3 columns: RA (h:m:s), Dec (d:m:s), S/N
    -p    PSF of a CB in fits format
    --o    Fractional sensitivity level at which CBs are tiled to overlap
    --r    Resolution of PSF in units of arcseconds per pixel
    --n    Number of beams to consider when creating overlap contours. Will
        pick the specified number of beams with the highest S/N values.
    --nsig Sets the number of standard deviation contours drawn.
    --s Draws known coordinates onto the plot for comparison.
    --scalebar Sets the length of the scale bar on the plot in arcseconds.
    --ticks Sets the spacing of ticks on the localisation plot.
    --clip Sets level below which CB PSF is set equal to zero.
    --zoom Automatically zooms in on the TABs.
    '''

    parser.add_argument('-f', dest='files', 
                nargs = '+', 
                type = str, 
                help="Detections files",
                required=True)
    parser.add_argument('-p',dest='psfs',
                nargs='+',
                type=str,
                help="PSF file",
                required=True)
    parser.add_argument('-c', dest='config', 
                nargs = 1, 
                type = str, 
                help="Configuration (json) file",
                required=False)
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
                help='Number of beams to use',
                default = [1000000])
    parser.add_argument('--nsig',dest='nsig',
                nargs = 1,
                type = int,
                help='Draws uncertainty contours up to this number of standard deviations.',
                default = [2])
    parser.add_argument('--s', dest='source',
                nargs = 1,
                type=str,
                help="Draws given coordinate location (format: hms,dms) on localisation plot",
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
                        default = [100],
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
    parser.add_argument('--fits', dest='fitsOut',
                        help = "Outputs .fits file of localisation region",
                        action = 'store_true')   

    options= parser.parse_args()

    return options

def readCoordsME(options):
    """
    Checks that arguments make sense, and if so returns data read from input file. 
    """
    RAs = []
    Decs = []
    
    data = np.concatenate([np.loadtxt(f,
                delimiter=' ',
                dtype=str,
                encoding="ascii") for f in options.files])

    # Sorting beams by S/N, so that Gaussian assumption holds better. 
    # You generally want LOW/HIGH S/N beam pairs.
    data = data[np.argsort(data[:,2])[::-1]] 

    c = SkyCoord(data[:,0], data[:,1], frame="icrs", unit=(u.hourangle, u.deg))
    best_cand = np.argsort(data[:,2])[-1]

    # Calculate number of beam pairs
    n_comb = math.factorial(len(data))/(math.factorial(2)*
                            math.factorial(len(data)-2)) 

    if options.source:
            [ra,dec] = options.source[0].split(',')
            b = SkyCoord(ra, dec, frame='icrs',unit=(u.hourangle, u.deg))
            boresightCoord = [b.ra.deg,b.dec.deg]
    else:
            bs_c = SkyCoord(data[:,0][best_cand],data[:,1][best_cand], frame='icrs', unit=(u.hourangle, u.deg))
            boresightCoord = [bs_c.ra.deg,bs_c.dec.deg]

    if options.overlap > 1.0 or options.overlap < 0:
            print("The OVERLAP parameter must be between 0 and 1")
            exit()

    elif options.npairs[0] > n_comb:
            options.npairs[0] = n_comb
            return data, c, boresightCoord
    else:
            return data, c, boresightCoord


def deg2pixME(w, c, psf, boresight, res):
    """
    Converts coordinates from degrees to pixels
    c must be a SkyCoord object and psf a numpy array.
    """

    coordsDeg = []
    for i in range(0, len(c)):
        coordsDeg.append([c.ra.deg[i], c.dec.deg[i]])

    ### Convert deg -> pix
    px = w.all_world2pix(coordsDeg, 1)
    c.ra.px = px[:, 0]
    c.dec.px = px[:, 1]

    return c


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    options = parseOptions(parser)

    f, ax = plt.subplots(figsize=(10,10))
    
    psf_ar = ut.readPSF(options.psfs[0], options.clipping[0])
    data, c, boresight = readCoordsME(options)
    c, w, array_width, array_height = co.deg2pix(c, psf_ar, boresight, options.res[0])

    for i in enumerate(options.files):
        options.file = [i[1]]
        print('\n' + options.file[0])
        options.psf = [options.psfs[i[0]]]
        data, c1, _ = ut.readCoords(options)

        psf_ar = ut.readPSF(options.psf[0], options.clipping[0])
        
        c1 = deg2pixME(w, c1, psf_ar, boresight, options.res[0])

        loglikelihood = np.zeros((array_height,array_width))
        loglikelihood += SK.make_map(array_height, array_width,
                             c1, psf_ar, options, data)

    if options.source:
        Splot.plot_known(w, options.source[0])
    
    Splot.make_ticks(ax, array_width, array_height, 
                    w, fineness=options.tickspacing[0])

    Splot.likelihoodPlot(f, ax, w, loglikelihood, options)

    if options.autozoom == True:
        ut.autozoom(ax, c, options)

    if options.fitsOut == True:
        likelihood = ut.norm_likelihood(loglikelihood)
        ut.write2fits(w, likelihood)

    plt.savefig(options.file[0]+'.png',dpi=300)
    plt.show()
