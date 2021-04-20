# SeeKAT multibeam localiser
Requires python 2

• Reads in list of detections (RA,Dec,S/N) and PSF. 

• Arranges beam PSFs appropriately relative to each other.

• Calculates localisation contours (+ 1 sigma errors) where ratio of PSFs match ratio of detections' S/N.

• Adds probability distributions for each pair of beams to create a localisation likelihood map.

Usage: python SeeKAT.py -f {coordinates file} -p {.fits file} --r {PSF resolution} --o {fractional overlap}

OR

python SeeKAT.py -f {coordinates file} -p {.fits file} --c {.json file} --r {PSF resolution}
