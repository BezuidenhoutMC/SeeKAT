# SeeKAT multibeam localiser
Requires python 2

• Reads in list of detections (RA,Dec,S/N) and PSF. 

• Arranges beam PSFs appropriately relative to each other.

• Calculates localisation contours (+ 1 sigma errors) where ratio of PSFs match ratio of detections' S/N.

• Adds probability distributions for each pair of beams to create a localisation likelihood map.

Usage: python SeeKAT.py -f {coordinates file} -p {.fits file} --r {PSF resolution} --o {fractional overlap}

OR

python SeeKAT.py -f {coordinates file} -p {.fits file} --c {.json file}

For more info see https://docs.google.com/document/d/102t90nL9pWn0I7aAAHnaN25Z9KS-6o4EVq6FerdN1Zc/edit
