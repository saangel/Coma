import glob as g
import subprocess as s
from astropy.io import fits as p
from joblib import Parallel,delayed
import os
import numpy as np
home=os.getcwd()
workdir="/Users/simon/Documents/Coma"
sextractor="sex"
det_filt="/usr/local/scisoft/packages/sextractor-2.5.0/config/gauss_1.5_3x3.conv"
os.chdir(workdir)
fits=g.glob('Data/final_drz_sci_imgs/*sci.fits')+g.glob('treasury/*sci*fits')
weights=g.glob('Data/final_drz_wht_imgs/*wht*fits')+g.glob('treasury/*wht*fits')
defparam1=open('default_1.param','w')
params1='NUMBER\nALPHA_J2000\nDELTA_J2000\nXWIN_IMAGE\nYWIN_IMAGE\nFLUX_RADIUS\nKRON_RADIUS\
\nMAG_APER\nMAG_AUTO\nELONGATION\nTHETA_IMAGE\nA_IMAGE\nB_IMAGE\nFLAGS\nFLUX_APER\n\
FLUXERR_APER\nSNR_WIN\nVIGNET(10,10)'
defparam1.write(params1)
defparam1.close()
defparam2=open('default_2.param','w')
params2='ALPHA_J2000\nDELTA_J2000\nFLUX_RADIUS\nFWHM_IMAGE\nELLIPTICITY\nFWHMPSF_IMAGE\
\nMAG_AUTO\nMAGERR_AUTO\nMAG_PSF\nMAGERR_PSF\nMAG_APER(5)\nMAGERR_APER(5)\nSPREAD_MODEL\
\nSPREADERR_MODEL\nFLUX_MAX\nA_IMAGE\nB_IMAGE\nNUMBER'
defparam2.write(params2)
defparam2.close()
zeros= {'F475W':26.0677371807,'F814W': 25.9365790112 }
#raw_input("break here")
def magic(f,w):
	hdu=p.open(f)
	'''
	try:
		w=weights[i]

	except:
		p.writeto(f[:-5]+'_sci.fits',hdu['SCI'].data,hdu['SCI'].header,clobber=True,output_verify='ignore')
		p.writeto(f[:-5]+'_wht.fits',hdu['WHT'].data,hdu['WHT'].header,clobber=True,output_verify='ignore')	
		hdu.close()
		f=f[:-5]+'_sci.fits'
		w=f[:-5]+'_wht.fits'
		hdu=p.open(f)
	'''	
	head=hdu[0].header
	filt=head['PHOTMODE'].split()[2]
	zp=-2.5*np.log10(float(head['PHOTFLAM']))+float(head['PHOTZPT'])-5*np.log10(float(head['PHOTPLAM']))+18.6921
	hdu.close()
	#return [filt, zp]
	if ~np.isfinite(zp):
		zp=zeros[filt]
	#print f,w,zp,filt,head['PHOTFLAM'],head['PHOTPLAM']
	
	cat='Catalogues/'+f.replace('/','_')[:-5].replace('.','')+'_'+filt+'.dat'
	print cat
	s.call([sextractor,f,'-PARAMETERS_NAME','default_1.param','-FILTER_NAME',det_filt,\
            '-CATALOG_NAME',cat,'-CATALOG_TYPE','FITS_LDAC','-VERBOSE_TYPE','QUIET',\
            "-DETECT_THRESH","1.5","-ANALYSIS_THRESH","1.5"])

	psf_model='PSFS/'+f.replace('/','_')[:-5].replace('.','')+'_'+filt+'.psf'
	s.call(['psfex',cat,'-PHOTFLUX_KEY','FLUX_APER','-PHOTFLUXERR_KEY','FLUXERR_APER'\
            ,'-PSFVAR_NSNAP','4','-PSF_SIZE','15,15','-CHECKPLOT_TYPE','NONE','-CHECKIMAGE_TYPE','NONE'\
            ,'-NTHREADS','1',"-CENTER_KEYS","XWIN_IMAGE,YWIN_IMAGE","-PSFVAR_KEYS"\
            ,"XWIN_IMAGE,YWIN_IMAGE",'-PSF_DIR','PSFS'])
	s.call(['rm',cat])

	cat_psf='Catalogues/'+f.replace('/','_')[:-5].replace('.','')+'_'+filt+'.dat'
	s.call([sextractor,f,'-PARAMETERS_NAME','default_2.param','-FILTER_NAME',det_filt,\
            '-CATALOG_NAME',cat_psf,'-CATALOG_TYPE','FITS_LDAC','-VERBOSE_TYPE','QUIET',\
            '-MAG_ZEROPOINT',str(zp),'-WEIGHT_TYPE','MAP_WEIGHT','-WEIGHT_IMAGE',w,\
            "-DETECT_THRESH","1.5","-ANALYSIS_THRESH","1.5","-PSF_NAME",psf_model,"PHOT_APERTURES","4,6,8,10,12"])
#	raw_input("Ready?")

Parallel(n_jobs=-1,verbose=10)(delayed(magic)(ff,ww) for ff,ww in zip(fits,weights))

# filts = ['F160W','F336W','F814W']
# for f in filts:
# 	fits=g.glob('*drz*'+f+'.fits')+g.glob('MAST_2016-06-13T2030/HST/*/*'+f+'.fits')+g.glob('../First_Batch/final_drz_sci_imgs/*160*'+f+'.fits')+g.glob('../First_Batch/final_drz_sci_imgs/*336*'+f+'.fits')+g.glob('../First_Batch/final_drz_sci_imgs/*814*'+f+'.fits')
# 	command = ['swarp']
# 	for fit in fits:
# 		command.append(fit)
# 	command.append('-IMAGEOUT_NAME')
# 	command.append('coadd_'+f+'.fits')
# 	s.call(command)


	
