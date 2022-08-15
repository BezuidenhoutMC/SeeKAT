# SeeKAT multibeam localiser

• Reads in list of detections (RA, Dec, S/N) and the beam PSF. PSFs can be modelled using MOSAIC (https://github.com/wchenastro/Mosaic).

• Computes a covariance matrix of the S/N value ratios, assuming 1-sigma Gaussian errors on each measurement.

• Models the aggregate beam response by arranging beam PSFs appropriately relative to each other.

• Calculates a likelihood distribution of obtaining the observed S/N in each beam according to the modelled response.

• Plots the likelihood function over RA and Dec with 1-sigma uncertainty, overlaid on the beam coordinates and sizes.

Usage: python SeeKAT.py -f {coordinates file} -p {.fits file} --r {PSF resolution} --o {fractional overlap}

OR

python SeeKAT.py -f {coordinates file} -p {.fits file} --c {.json file} --r {PSF resolution}

Customisation options:

--n Computes the likelihood map using only the n brightest pairs of beams.

--clip All values of the CB PSF below this value are set to zero. Helps negate low-level sidelobes.

--s In the format RA(hms),Dec(dms) adds a marker for known coordinates to plot.

--scalebar Sets the length of the scalebar on the plot in arcseconds. Set to 0 to omit it altogether.

--zoom Automatically zooms in on the TABs.

--ticks Sets the spacing between axis ticks in number of pixels.

--nsig Draws uncertainty contours up to this number of standard deviations.

--fits Writes likelihood to fits file.
