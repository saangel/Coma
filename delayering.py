import numpy as np
import glob as g
import aplpy
from astropy.io import fits as p
import os

home=os.getcwd()
os.chdir("/Users/simon/Documents/Coma")
fits=g.glob('treasury/j*drc*fits')+g.glob('treasury/j*drz*fits')
filters=["F336W","F475W","F814W","F160W"]
for f in fits:
  hdu=p.open(f)
  if len(hdu) > 1:
    try:
      filt = hdu[0].header["PHOTMODE"].split()[2]
      head=hdu[0].header
    except:
      filt = hdu["SCI"].header["PHOTMODE"].split()[2]
      head=hdu["SCI"].header
    if filt in filters:
        p.writeto(f[:-5]+"_"+filt+'_sci.fits',hdu['SCI'].data,head,clobber=True,output_verify='ignore')
        p.writeto(f[:-5]+"_"+filt+'_wht.fits',hdu['WHT'].data,head,clobber=True,output_verify='ignore')
        hdu.close()
  hdu.close()

fits_del = g.glob("*drc_sci.fits")+g.glob("*drz_sci.fits")+g.glob("*drc_wht.fits")+g.glob("*drz_wht.fits")
for f in fits_del:
  os.remove(f)