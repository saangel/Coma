import glob as g
import subprocess as s
from astropy.io import fits as p
from joblib import Parallel,delayed
from astropy.wcs import WCS
import numpy as np
import os, platform, copy

home=os.getcwd()
workdir="/home/simon/Documents/Coma/First_Batch/final_drz_sci_imgs/"
sextractor="sextractor"
filt='/usr/share/sextractor/gauss_1.5_3x3.conv'
if platform.system() == "Darwin":
	workdir="/Users/simon/Documents/Coma/First_Batch/final_drz_sci_imgs/"
	sextractor="sex"
	filt="/usr/local/scisoft/packages/sextractor-2.5.0/config/gauss_1.5_3x3.conv"
os.chdir(workdir)
params='NUMBER\nALPHA_J2000\nDELTA_J2000\nXWIN_IMAGE\nYWIN_IMAGE\nFLUX_RADIUS\nKRON_RADIUS\nMAG_APER\nMAG_AUTO\nELONGATION\nTHETA_IMAGE\nA_IMAGE\nB_IMAGE\nFLAGS\nFLUX_AUTO'
defparam=open('default_galfit.param','w')
defparam.write(params)
defparam.close()
files=g.glob("*160*sci.fits")#+g.glob("*336*sci.fits")+g.glob("*475*sci.fits")+g.glob("*814*sci.fits")
skips=1
for f in files[::skips]:
	cat_name=f[:-5]+".cat"
	checkimage_name=f[:-5]+"_seg.fits"
	hdu=p.open(f)
	head=hdu[0].header
	pix_scale=head["D001SCAL"]
	zp=-2.5*np.log10(float(head['PHOTFLAM']))+float(head['PHOTZPT'])-5*np.log10(float(head['PHOTPLAM']))+18.6921
	im_data=hdu[0].data
	im_size=hdu[0].header['NAXIS1'],hdu[0].header['NAXIS2']
	s.call([sextractor,f,'-PARAMETERS_NAME','default_galfit.param','-FILTER_NAME',filt,'-CATALOG_NAME',cat_name,'-CATALOG_TYPE','FITS_LDAC','-VERBOSE_TYPE','QUIET','-MAG_ZEROPOINT',str(zp),'-CHECKIMAGE_TYPE','SEGMENTATION','-CHECKIMAGE_NAME',checkimage_name])
	cat_data=p.open(cat_name)[2].data
	a_image=cat_data["A_IMAGE"]
	b_image=cat_data["B_IMAGE"]
	border=min(im_size)
	a_limit=0.25
	if "160" in f: a_limit=1
  #  print a_limit*pix_scale
	galaxy_condition=((a_image*pix_scale>=a_limit)&(b_image*pix_scale>=a_limit))&(cat_data["MAG_AUTO"]<99)&(cat_data["B_IMAGE"]/cat_data["A_IMAGE"]>0.1)&(cat_data["XWIN_IMAGE"]<.99*border)&(cat_data["YWIN_IMAGE"]<.99*border)
	galaxy_like=cat_data[galaxy_condition]
	not_gal_like=cat_data[~galaxy_condition]
	gal_nums=galaxy_like["NUMBER"]
	not_gal_nums=not_gal_like["NUMBER"]
  #  print pix_scale
	print f,pix_scale
  #  print len(galaxy_like)
  
	############################################
	#
	#        Reverse masking process -- done!!
	#
	############################################
	gal_feedme_single="galfit_"+f[:-5]+"_single.feedme"
	gal_feedme_double="galfit_"+f[:-5]+"_double.feedme"
	seg_image=f[:-5]+"_seg_image.fits"
	'''
	gal_file_single=open(gal_feedme_single,"w")
	gal_file_double=open(gal_feedme_double,"w")
	gal_file_single.write("A) "+seg_image+"\n")
	gal_file_single.write("B) "+f[:-5]+"_gal_single.fits\n")
	gal_file_single.write("C) none\n")
	gal_file_single.write("H) 1 %.0f 1 %.0f \n" % (im_size[0],im_size[1]))
	gal_file_single.write("I) %.0f %.0f\n" % (np.round(im_size[0]/5,-2),np.round(im_size[1]/5,-2)))
	gal_file_single.write("J) %.6f 1 \n" %zp)
	gal_file_single.write("K) %.6f %.6f 1 1 \n" % (pix_scale,pix_scale))
	gal_file_single.write("O) regular\n")
	gal_file_single.write("P) 0\n")
	gal_file_double.write("A) "+seg_image+"\n")
	gal_file_double.write("B) "+f[:-5]+"_gal_double.fits\n")
	gal_file_double.write("C) none\n")
	gal_file_double.write("H) 1 %.0f 1 %.0f \n" % (im_size[0],im_size[1]))
	gal_file_double.write("I) %.0f %.0f\n" % (np.round(im_size[0]/5,-2),np.round(im_size[1]/5,-2)))
	gal_file_double.write("J) %.6f 1 \n" %zp)
	gal_file_double.write("K) %.6f %.6f 1 1 \n" % (pix_scale,pix_scale))
	gal_file_double.write("O) regular\n")
	gal_file_double.write("P) 0\n")

	gal_file_single.close()
	gal_file_double.close()
	
	'''
	#if head['BUNIT'] == "ELECTRONS/S": 
	#	im_data=im_data*head['EXPTIME']
	#	head['BUNIT'] == "ELECTRONS"
	
	gal_file_name="galfit_"+f[:-5]+".feedme"
	gal_file=open(gal_file_name,"w")
	gal_file.write("A) "+seg_image+"\n")
	gal_file.write("B) "+f[:-5]+"_gal.fits\n")
	gal_file.write("C) none\n")
	gal_file.write("H) 1 %.0f 1 %.0f \n" % (im_size[1],im_size[0]))
	gal_file.write("I) %.0f %.0f\n" % (np.round(im_size[1]/20,-2),np.round(im_size[0]/20,-2)))
	gal_file.write("J) %.6f 1 \n" %zp)
	gal_file.write("K) %.6f %.6f 1 1 \n" % (pix_scale,pix_scale))
	gal_file.write("O) regular\n")
	gal_file.write("P) 0\n")
	file_list=[]
	for j,g in enumerate(galaxy_like):
		
		gal_file.write("0) sersic\n")
		gal_file.write("1) %.2f %.2f 1 1\n" % (g["XWIN_IMAGE"],g["YWIN_IMAGE"]))
		gal_file.write("3) %.2f  1 \n" % (zp-2.5*np.log10(g["FLUX_AUTO"])))
		gal_file.write("4) %.2f  1 \n" % g["FLUX_RADIUS"])
		gal_file.write("5) 4 1\n")
		gal_file.write("9) %.3f 1\n" % (g["B_IMAGE"]/g["A_IMAGE"]))
		gal_file.write("10) %.2f 1\n" % (g["THETA_IMAGE"]+90))
		gal_file.write("\n")
		'''
		gal_file.write("0) sersic\n")
		gal_file.write("1) %.2f %.2f 1 1\n" % (g["XWIN_IMAGE"],g["YWIN_IMAGE"]))
		gal_file.write("3) %.2f  1 \n" % (zp-2.5*np.log10(3*g["FLUX_AUTO"]/4)))
		gal_file.write("4) %.2f  1 \n" % g["FLUX_RADIUS"])
		gal_file.write("5) 4 1\n")
		gal_file.write("9) %.3f 1\n" % (g["B_IMAGE"]/g["A_IMAGE"]))
		gal_file.write("10) %.2f 1\n" % (g["THETA_IMAGE"]+90))	
		'''

	gal_file.close()	
	'''
	
	data_seg=p.open(checkimage_name)[0].data
	mastermask=np.zeros_like(data_seg,dtype=bool)
	im_data_2=copy.deepcopy(im_data)
	for n in not_gal_nums:
		mask = data_seg==n 
		im_data_2=np.where(~mask,im_data_2,0)
		print float(n)/len(not_gal_nums)

	p.writeto(seg_image,im_data_2,head,clobber=True,output_verify='ignore')
	'''
	hdu.close()
	
	#for fl in file_list:
	#	s.call("/Users/simon/Documents/Galfit/galfit "+fl,shell=True)	
	#s.call("/Users/simon/Documents/Galfit/galfit "+gal_file_name,shell=True)
	#s.call("/Users/simon/Documents/Galfit/galfit "+gal_feedme_double,shell=True)
	#plt.matshow(np.where(mastermask,im_data,0),norm=LogNorm(vmin=0.001, vmax=1))
	#plt.show()
	#plt.matshow(im_data-np.where(mastermask,im_data,0),norm=LogNorm(vmin=0.001, vmax=1))
	#plt.show()
	#plt.matshow(im_data-np.where(~mastermask,im_data,0),norm=LogNorm(vmin=0.001, vmax=1))
	#plt.show()
	
	#s.call("/Users/simon/Documents/Galfit/galfit "+gal_feedme_single+' -imax 200',shell=True)








os.chdir(home)
