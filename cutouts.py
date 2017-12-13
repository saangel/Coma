import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.colors import LogNorm
from sklearn.neighbors import KernelDensity


filters=['F160W']#,'F336W']
vmax=[25,0.6]


test=fits.open('../Catalogues/Compact_objects_unfiltered_3008.fits')[-1]
data=test.data#np.concatenate((test.data[test.data['TC_rowSubsetFlags']==2],test.data[test.data['TC_rowSubsetFlags']==3]))
nobj=len(data)
grid=int(np.ceil(np.sqrt(nobj))),int(np.ceil(np.sqrt(nobj)))
end=False
ii=0
size=15

#white not end:
for i,f in enumerate(filters):
	ra_cat,dec_cat=data['ALPHA_J2000_'+str(i+2)],data['DELTA_J2000_'+str(i+2)]
	image=fits.open('../FITS files/coadd_'+f+'.fits')[0]
	imdata=image.data

	wcs = WCS(image.header)
	y,x=wcs.all_world2pix(ra_cat,dec_cat,1)
#	cutouts_data=np.array([imdata[x[i]-size:x[i]+size,y[i]-size:y[i]+size] for i in range(len(x))])
#	plt.hist(cutouts_data.flatten(),bins=100)
#	plt.yscale('log')
#	plt.show()
#	vmax=float(raw_input('Ingrese el maximo para graficar'))
	z = 0
	for j in range(grid[0]):
		for k in range(grid[1]):
			ax=plt.subplot2grid((grid[0],grid[1]+1),(j,k))
			ax.axis('off')
			if z<nobj:
				fig=ax.matshow(imdata[x[z]-size:x[z]+size,y[z]-size:y[z]+size],cmap='jet',norm=LogNorm(vmin=0.001, vmax=vmax[i]))
				z=z+1
		cax=plt.subplot2grid((grid[0],grid[1]+1),(0,k+1),rowspan=grid[0])
	cb=plt.colorbar(fig,cax=cax)
	cb.set_label("Counts")
	plt.savefig('test_.pdf')
	plt.clf()
	#image.close()

#for i in range(len(x)):
	#plt.imshow(imdata[x[i]-size:x[i]+size,y[i]-size:y[i]+size],cmap='jet',interpolation='bicubic')
	#plt.colorbar()
	#plt.show()
