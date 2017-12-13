'''Script to generate coverage maps and pointing coordinates
for an observing plan with the VIMOS imaging camera.'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.transforms as t
import glob as g
from astropy.io import fits as p
import sys

def draw_chips(x,y,chip_x0,chip_y,gap_sht0,angle,colour,fil):
	'''Generate rectangular patches to display on the
	coverage map.'''
	y0 = y-chip_y/2
	chip_x=chip_x0
	gap_sht = 0

	x0 = x-chip_x/2

	tr = t.Affine2D().rotate_deg_around(x,y,angle)
	tf=tr+ts
	if fil==1:
		rect = patches.Rectangle([x0,y0],chip_x,chip_y,edgecolor="none",facecolor=colour,alpha=.25)
	else:
		rect = patches.Rectangle([x0,y0],chip_x,chip_y,edgecolor=colour,facecolor="none",alpha=.15,ls="-.")
	rect.set_transform(tf)

	return rect

plt.clf()
fig, (ax) = plt.subplots(1,1,figsize=(6,6))
ts=ax.transData

xlimmin = []
xlimmax = []
ylimmin = []
ylimmax = []
filters=["F336W","F475W","F814W","F160W"]
colours=["blue","green","yellow","red"]
scale = {'F475W':0.05,'F814W': 0.05, 'F336W': 0.04, 'F160W': 0.13}

posx=[]
posy=[]
for f,c in zip(filters,colours):
	#visits = g.glob("/Users/simon/Documents/Coma/Data/*sci*/*"+f+"*sci*") #+ g.glob("/Users/simon/Documents/Coma/treasury/*"+f+"*sci.fits")
	visits_archive = g.glob("/Users/simon/Documents/Coma/treasury/*"+f+"*sci.fits")
	visits_new = g.glob("/Users/simon/Documents/Coma/Data/*sci*/*"+f+"*sci*")
	visits=visits_archive+visits_new
	filled=np.concatenate((np.zeros_like(visits_archive,dtype=np.int),np.ones_like(visits_new,dtype=np.int)))
	#print len(visits)
	for j,v in enumerate(visits):
		hdu=p.open(v)
		head=hdu["PRIMARY"].header
		centre_ra = head["crval1"]
		centre_dec = head["crval2"]

		try:
			chip_x = head["NAXIS1"] #size of the chip along the x axis in pixels
		except:
			chip_x = head["SIZAXIS1"]
		try: 
			chip_y = head["NAXIS2"] #size of the chip along the y axis in pixels
		except:
			chip_y = head["SIZAXIS2"]
		posx.append(centre_ra)
		posy.append(centre_dec)
		try:
			pres = head["D001SCAL"]
		except:
			pres = scale[f]
		angle=head["orientat"]
		gap_sht = 0. #degrees
		gap_lng = 0. #degrees
		chip_x = chip_x*pres/3600. #degrees
		chip_y = chip_y*pres/3600. #degrees
		xlimmin.append(centre_ra-chip_x/np.cos(centre_dec*np.pi/180.))
		xlimmax.append(centre_ra+chip_x/np.cos(centre_dec*np.pi/180.))
		ylimmin.append(centre_dec-chip_y)
		ylimmax.append(centre_dec+chip_y)

		ax.add_patch(draw_chips(centre_ra,centre_dec,chip_x,chip_y,gap_sht,-angle,c,filled[j]))
		print pres    
#	print filled
	print "plotting "+f 


ax.set_xlim(195.2,194.8)
ax.set_ylim(27.8,28.2)
ax.set_xlabel(r'RA [deg]')
ax.set_ylabel(r'Dec [deg]')
plt.scatter(194.8987875,27.9593889,marker='^',color='black',s=50,label='NGC4874')
plt.scatter(195.0337375,27.9770250,marker='v',color='black',s=50,label='NGC4889')
plt.legend(loc='best',ncol=2)
#plt.savefig('deleteme_dithering_pattern.pdf',bbox_inches='tight')

# plt.savefig('dithering_pattern.pdf')

plt.show()
plt.clf()
	