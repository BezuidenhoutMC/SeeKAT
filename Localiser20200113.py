#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from astropy.coordinates import SkyCoord
import astropy.units  as u
import matplotlib.patches as ptch
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', dest='file', nargs = 1, type = str, help="Candidate file",required=True)
options= parser.parse_args()


pulse = np.genfromtxt(options.file[0],delimiter='\t',dtype=None,names=["files","MJD","DM","SN","beam","RA","Dec","fil"],encoding="ascii")

min_ra=280.85
max_ra=280.92
min_dec=-7.995
max_dec=-7.925


#----------
#handling files with one candidate
try:
	if pulse["MJD"][0] < 58794.63333333:
		overlap=0.95
	else:
		overlap=0.98

except:
	if float(pulse["MJD"]) < 58794.63333333:
		overlap=0.95
		el_width = 0.00156
		el_height=0.00166
		el_angle=-4.553
		psf='psf95.fits'
	else:
		overlap=0.98
		el_width = 0.00248
		el_height=0.00265
		el_angle=-4.553
		psf='psf98.fits'

	ax = plt.gca()
	co = SkyCoord(str(pulse["RA"]),str(pulse["Dec"]), frame="icrs",unit=(u.hourangle,u.deg))
	sc=plt.scatter(co.ra.deg,co.dec.deg,zorder=4)
	plt.xlim(min_ra,max_ra)
	plt.ylim(min_dec,max_dec)
	e=ptch.Ellipse(xy=(co.ra.deg,co.dec.deg),width=2*el_width,height=2*el_height,angle=el_angle,fill=False,color='white')
	ax.add_artist(e)
	plt.grid(linestyle='-')
	plt.title(options.file[0])
	filename_split = str(pulse["files"]).split('/')
	plt.xlabel("Best candidate: S/N="+str(pulse["SN"])+" in node "+filename_split[0][6:8]+ "\n" + filename_split[5]+'/'+filename_split[6]+'/'+str(pulse["fil"]))
	plt.tight_layout()
	plt.savefig("Contours_"+str(options.file[0])+".png")
	exit()

#-----------------

if overlap==0.98:
	el_width = 0.00156
	el_height=0.00166
	el_angle=-4.553
	psf='psf95.fits'
if overlap==0.95:
	el_width = 0.00248
	el_height=0.00265
	el_angle=-4.553
	psf='psf98.fits'

hdul = fits.open(psf)
hdr = hdul[0].header
data = hdul[0].data



#--------- Plot pulse coords ---------

co = SkyCoord(pulse["RA"],pulse["Dec"], frame="icrs",unit=(u.hourangle,u.deg))
extent=(min(co.ra.deg),max(co.ra.deg),min(co.dec.deg),max(co.dec.deg))
best_cand = np.argsort(pulse["SN"])[-1]
best_coords=(co.ra.deg[best_cand],co.dec.deg[best_cand])
box_c = [280.8874,-7.9601]

c=co

c.ra.deg *= (63/(2*el_width))
min_ra *= (63/(2*el_width))
max_ra *= (63/(2*el_width))
box_c[0] *= (63/(2*el_width))
min_ra = min_ra - c.ra.deg[best_cand] + 525
max_ra = max_ra - c.ra.deg[best_cand] + 525
box_c[0] -= c.ra.deg[best_cand]-525
c.ra.deg -= c.ra.deg[best_cand]-525

c.dec.deg *= (63/(2*el_height))
box_c[1] *= (63/(2*el_height))
min_dec *= (63/(2*el_height))
max_dec *= (63/(2*el_height))
min_dec -= c.dec.deg[best_cand]-525
max_dec -= c.dec.deg[best_cand]-525
box_c[1] -= c.dec.deg[best_cand]-525
c.dec.deg -= c.dec.deg[best_cand]-525

sc=plt.scatter(c.ra.deg,c.dec.deg,c=pulse["SN"],zorder=4)
cb = plt.colorbar(sc)
cb.set_label('S/N')
#plt.scatter(np.median(c.ra.deg),np.median(c.dec.deg),marker='x',color='gray',zorder=6)
ax = plt.gca()

#for i in range(0,len(c.ra.deg)):
#	e=ptch.Ellipse(xy=(c.ra.deg[i],c.dec.deg[i]),width=63,height=63,angle=el_angle,fill=False,color='white')
#	ax.add_artist(e)

#-----BOX--------
#plt.scatter(box_c[0],box_c[1],marker='x',color='magenta',zorder=9999)
width = 0.002
height = 0.002
width *= (63/(2*el_width))
height *= (63/(2*el_height))
box_c[0]-=0.5*width
box_c[1]-=0.5*height
rect=ptch.Rectangle(xy=box_c,width=width,height=height,linewidth=1,edgecolor='magenta',facecolor='none',zorder=9999)
#ax.add_artist(rect)
#---------------


#-----------Interpolation--------------------

data*=max(pulse["SN"])
points = np.zeros((400,2))
values = []
n=0
for i in range (0,20):
	for j in range(0,20):
		points[n] = [i, j]
		values.append(data[i,j])
		n+=1


grid_x,grid_y = np.mgrid[0:18.98:1000j,0:18.98:1000j]
model_prime = griddata(points,values,(grid_x,grid_y), method='cubic')

#------------Ratio Contours--------------
model_prime[model_prime<0.05*max(pulse["SN"])] = np.nan

n,m = model_prime.shape
model_padded = np. zeros((9*n,9*m))			#beam8
model_padded[4*n:5*n,4*m:5*m]=model_prime


plt.gca()
psf = plt.imshow(model_prime,vmin=8)


for i in range (0,len(pulse)):
	model = model_padded[4*n-int(c.dec.deg[i]-c.dec.deg[best_cand]):5*n-int(c.dec.deg[i]-c.dec.deg[best_cand]),4*m-int(c.ra.deg[i]-c.ra.deg[best_cand]):5*m-int(c.ra.deg[i]-c.ra.deg[best_cand])]
	plt.contour(model,colors='white',levels=[0.98*max(pulse["SN"])])#,extent=extent)
	
	if len(pulse)>6:	
		if pulse["SN"][i]>=pulse["SN"][np.argsort(pulse["SN"])[-6]]:
			error = (1/max(pulse["SN"]) + 1/pulse["SN"][i])
			psf = plt.contour(np.divide(model_prime,model),colors='red',levels=[max(pulse["SN"])/pulse["SN"][i]],zorder=800)#,extent=extent)
			#ERROR CONTOURS
			plt.contour(np.divide(model_prime,model),linestyles='dashed',linewidths=0.7,colors='magenta',levels=[max(pulse["SN"])/pulse["SN"][i]] - error)
			plt.contour(np.divide(model_prime,model),linestyles='dashed',linewidths=0.7,colors='magenta',levels=[max(pulse["SN"])/pulse["SN"][i]] + error)
	else:
		error = (1/max(pulse["SN"]) + 1/pulse["SN"][i])
		psf = plt.contour(np.divide(model_prime,model),colors='red',levels=[max(pulse["SN"])/pulse["SN"][i]])#,extent=extent)
		#ERROR CONTOURS
		plt.contour(np.divide(model_prime,model),linestyles='dashed',linewidths=0.7,colors='magenta',levels=[max(pulse["SN"])/pulse["SN"][i]] - error)
		plt.contour(np.divide(model_prime,model),linestyles='dashed',linewidths=0.7,colors='magenta',levels=[max(pulse["SN"])/pulse["SN"][i]] + error)
#plt.axis([200,900,200,900])

locs = np.arange(-300,1300,50)
labels = [round((float(item)-525)*(2*el_width/63)+best_coords[0],4) for item in locs]
plt.xticks(locs, labels,rotation=90)
#plt.vlines(x=525,ymin=0,ymax=1000,colors='magenta')
labels = [round((float(item)-525)*(2*el_height/63)+best_coords[1],5) for item in locs]
plt.yticks(locs, labels)
#plt.hlines(y=525,xmin=0,xmax=1000,colors='magenta')
plt.grid(linestyle='-')
plt.axis([min_ra,max_ra,min_dec,max_dec])
plt.title(options.file[0])



filename_split = str(pulse["files"][best_cand]).split('/')
plt.xlabel("Right Ascension (deg) \nBest candidate: S/N="+str(pulse["SN"][best_cand])+" in node "+filename_split[0][6:8]+ "\n" + filename_split[5]+'/'+filename_split[6]+'/'+str(pulse["fil"][best_cand]))
plt.ylabel("Declination (deg)")
plt.tight_layout()
plt.savefig("Contours_"+str(options.file[0])+".png",dpi=300)
plt.show()
plt.close()
#---------------------------------------
#RESIZING



