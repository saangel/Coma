from astropy.io import fits as p
from astropy.wcs import WCS
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from glob import glob as g
from copy import deepcopy as dc
import numpy as np
import os
import subprocess as s

home=os.getcwd()
workdir="/home/simon/Documents/Coma/Data"
sextractor="sextractor"
det_filt="/usr/share/sextractor/gauss_1.5_3x3.conv"
figdir="/home/simon/Documents/Coma/coding/figures/apercorr/"
napers=30
a0=1
sep=1
apers=range(a0,sep*napers+a0,sep)
apers_string=""
for aper in apers:
	apers_string+=str(aper)+","
apers_string=apers_string[:-1]
os.chdir(workdir)
parfilename="apercorr.param"
paramfile=open(parfilename,"w")
paramfile.write("ALPHA_J2000\nDELTA_J2000\nXWIN_IMAGE\nYWIN_IMAGE\nMAG_AUTO \
	\nFLUX_RADIUS\nMAG_APER("+str(len(apers))+")\nMAGERR_APER("+str(len(apers))+")\
	\nNUMBER\nFLAGS\nELONGATION\nFWHM_IMAGE")
paramfile.close()
files=sorted(g("v*_sci.fits"))
skips=1
for f in files[::skips]:
	catname=f[:10]+"apercorr.fits"
	#print " ".join([sextractor,f+"[0]","-PARAMETERS_NAME",parfilename,"-FILTER_NAME",det_filt,"-PHOT_APERTURES",apers_string,
	#"-CATALOG_TYPE","FITS_LDAC","-CATALOG_NAME",catname])
	#s.call([sextractor,f+"[0]","-PARAMETERS_NAME",parfilename,"-FILTER_NAME",det_filt,"-PHOT_APERTURES",apers_string,
	#"-CATALOG_TYPE","FITS_LDAC","-CATALOG_NAME",catname,"-BACKPHOTO_TYPE","LOCAL"])
	hdu=p.open(catname)
	data=hdu[2].data
	data=data[data["FLUX_RADIUS"]>1]
	data=data[data["MAG_AUTO"]<50]
	od=dc(data)
	data=data[data["FLAGS"]==0]
	data=data[data["FLUX_RADIUS"]<2]
	selected_data=np.sort(data,order="MAG_AUTO")[2:5]
	for sd in selected_data:
		plt.plot(apers,sd["MAG_APER"].reshape(napers)-sd["MAG_APER"][5],"-",
		label=str(sd["NUMBER"])+" "+str(sd["FLAGS"])+" "+str(sd["ELONGATION"])+" "+str(sd["FWHM_IMAGE"]))
	plt.legend(loc="best")
	plt.savefig(figdir+f[:10]+"apercorr.png",density=300)
	plt.clf()
	plt.scatter(od["FLUX_RADIUS"],od["MAG_AUTO"],marker="o",s=5,c="grey")
	plt.scatter(selected_data["FLUX_RADIUS"],selected_data["MAG_AUTO"],marker="*",s=8,c="red")
	plt.savefig(figdir+f[:10]+"sample.png",density=300)
	plt.clf()
	hdu.close()
	imhdu=p.open(f)
	imdata=imhdu[0].data
	wcs = WCS(imhdu[0].header)
	fig,(ax1,ax2,ax3)=plt.subplots(1,3,sharex=True,sharey=True)
	axes=[ax1,ax2,ax3]
	for i,sd in enumerate(selected_data):
		y,x=wcs.all_world2pix(sd["ALPHA_J2000"],sd["DELTA_J2000"],1)
		axes[i].axis("off")
		axes[i].matshow(imdata[int(x)-napers:int(x)+napers,int(y)-napers:int(y)+napers],cmap='viridis',norm=LogNorm(vmin=0.001,vmax=1))
		axes[i].set_title(" N "+str(sd["NUMBER"]))
	plt.savefig(figdir+f[:10]+"cutouts.png",density=300)
	plt.clf()
	imhdu.close()
		
	
os.chdir(home)


