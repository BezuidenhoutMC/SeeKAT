# SeeKAT multibeam localiser
Requires python 2

• Reads in list of detections (RA,Dec,S/N). 

• Reads in configuration file with boresight coordinates, time of observation, and antenna configuration. Then generates the MeerKAT PSF using mosaic.

• Arranges beam PSFs appropriately relative to each other.

• Calculates localisation contours (+ 1 sigma errors) where ratio of PSFs match ratio of detections' S/N.

Usage: python SeeKAT.py -f {data file} -c {config file} -options

Options: -p {psf file}: read in pre-rendered PSF instead of generating a new one

Example data file (.txt) and config file (.cfg) included.
