from astropy.io import fits as p
from acstools import calacs
from glob import glob as g
from os.path import isfile
import subprocess as s

jref="/raid/coma/jref/"
ljref="ftp://ftp.stsci.edu/cdbs/jref/"
jref_keys=['BPIXTAB','CCDTAB','OSCNTAB','BIASFILE','PCTETAB','FLSHFILE','CRREJTAB','SHADFILE','PCTETAB','DRKCFILE','DARKFILE','PFLTFILE','IDCTAB','DFLTFILE','LFLTFILE','PHOTTAB','DGEOFILE','MDRIZTAB','CFLTFILE','SPOTTAB','GRAPHTAB','COMPTAB','IMPHTTAB','D2IMFILE','NPOLFILE','SNKCFILE','MLINTAB']


rawfiles=g("*_raw.fits")
for r in rawfiles:
	hdu=p.open(r)
	print r
	head=hdu[0].header
	for j in jref_keys:
		f=head[j]
		if f == "N/A":
			print "no file for "+j
		else:
			f=f.split("$")[1]
			if not isfile(jref+f):
				print f," doesnt exist"
				s.call(["wget",ljref+f])
				s.call(["mv",f,jref+f])
			else:
				print f," exists"

