<!-- # SeeKAT multibeam localiser -->
If you use SeeKAT for a publication, please cite my paper: https://ui.adsabs.harvard.edu/abs/2023RASTI...2..114B

NB: A browser-based interactive app version of SeeKAT is now available here: https://github.com/BezuidenhoutMC/SeeKAT-app/

• Reads in list of detections (RA, Dec, S/N) and the beam PSF. PSFs can be modelled using MOSAIC (https://github.com/wchenastro/Mosaic).

• Computes a covariance matrix of the S/N value ratios, assuming 1-sigma Gaussian errors on each measurement.

• Models the aggregate beam response by arranging beam PSFs appropriately relative to each other.

• Calculates a likelihood distribution of obtaining the observed S/N in each beam according to the modelled response.

• Plots the likelihood function over RA and Dec with 1-sigma uncertainty, overlaid on the beam coordinates and sizes.

<img src="https://user-images.githubusercontent.com/22096485/184672571-49ff4929-5ccf-4940-bf4c-03feb5e6b163.png" width="400">

<ins>Usage</ins>: python SeeKAT.py -f {coordinates file} -p {.fits file} --r {PSF resolution} --o {fractional overlap}

OR

python SeeKAT.py -f {coordinates file} -p {.fits file} --c {.json file} --r {PSF resolution}

<ins>Customisation options</ins>:

<b>--n</b> Computes the likelihood map using only the n brightest pairs of beams.

<b>--clip</b> All values of the CB PSF below this value are set to zero. Helps negate low-level sidelobes.

<b>--s</b> In the format RA(hms),Dec(dms) adds a marker for known coordinates to plot.

<b>--scalebar</b> Sets the length of the scalebar on the plot in arcseconds. Set to 0 to omit it altogether.

<b>--zoom</b> Automatically zooms in on the TABs.

<b>--ticks</b> Sets the spacing between axis ticks in number of pixels.

<b>--nsig</b> Draws uncertainty contours up to this number of standard deviations.

<b>--fits</b> Writes likelihood to fits file.
