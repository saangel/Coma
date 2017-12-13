import numpy as np
import glob as g
from astropy.io import fits as p
import os

home=os.getcwd()
gen=os.walk(home)
fits_dirs=gen.next()[1]
filters=["F336W","F475W","F814W","F160W"]

for f in fits_dirs:
  os.chdir(f)
  fit=g.glob("*drc.fits")
  if not fit:
	  fit=g.glob("*drz.fits")
  if not fit:
	  os.chdir(home)
	  continue
  fit=fit[0]
  hdu=p.open(fit)
 #print hdu.info()
  try:
	  filt= hdu["SCI"].header["PHOTMODE"].split(" ")[2]
  except:
	  filt= hdu[0].header["PHOTMODE"].split(" ")[2]
  #    head=hdu["SCI"].header
  #  if filt in filters:
  try:
	  print hdu["sci"].data*1
	  p.writeto(fit[:9]+"_"+filt+'_sci.fits',hdu['SCI'].data,hdu['SCI'].header,overwrite=True,output_verify='ignore')
	  p.writeto(fit[:9]+"_"+filt+'_rms.wht.fits',1/np.sqrt(hdu['WHT'].data*0.77),hdu["WHT"].header,overwrite=True,output_verify='ignore')	
  except:
	  print fit, "fail"
  #if not hdu['SCI'].data.all():
	#print fit, "fail"
	#os.chdir(home)
	#continue
  #else:
	  
  hdu.close()
  
  os.chdir(home)
'''
fits_del = g.glob("*drc_sci.fits")+g.glob("*drz_sci.fits")+g.glob("*drc_wht.fits")+g.glob("*drz_wht.fits")
for f in fits_del:
  os.remove(f)
'''
