# SeeKAT multibeam localiser

• Reads in list of detections (RA,Dec,S/N) and PSF. 

• Arranges beam PSFs appropriately relative to each other.

• Calculates localisation contours (+ 1 sigma errors) where ratio of PSFs match ratio of detections' S/N.

• Adds probability distributions for each pair of beams to create a localisation likelihood map.

Usage: python SeeKAT.py -f {coordinates file} -p {.fits file} --r {PSF resolution} --o {fractional overlap}

OR

python SeeKAT.py -f {coordinates file} -p {.fits file} --c {.json file} --r {PSF resolution}

Customisation options:

--n Computes the likelihood map using only the n brightest pairs of beams.

--clip All values of the CB PSF below this value are set to zero. Helps negate low-level sidelobes.

--s In the format RA,Dec adds a marker for known coordinates to plot

--scalebar Sets the length of the scalebar on the plot in arcseconds. Set to 0 to omit it altogether.

--zoom Automatically zooms in on the TABs.

--ticks Sets the spacing between axis ticks in number of pixel.
