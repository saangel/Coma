from matplotlib.markers import MarkerStyle as mk
markers=mk.filled_markers
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import glob as g
import numpy as np
from astropy.coordinates import SkyCoord as SC
from astropy import units as u
from astropy.io import ascii
cats=sorted(g.glob('Catalogues/v*dat'))
visits=np.array(cats).reshape((len(cats)/4,4))
color=cm.jet(np.linspace(0,1,len(visits)))
markers=markers[:len(visits)]
for v,c in zip(visits,color):
	alpha,delta=np.loadtxt(v[0],usecols=(0,1),unpack=True)
	plt.scatter(alpha,delta,color=c,s=2,label=str(v[0]).split('/')[1][:3])
	alpha,delta=np.loadtxt(v[2],usecols=(0,1),unpack=True)
	plt.scatter(alpha,delta,color=c,s=2)
plt.scatter(194.8987875,27.9593889,marker='^',color='black',s=50,label='NGC4874')
plt.scatter(195.0337375,27.9770250,marker='v',color='black',s=50,label='NGC4889')
data=ascii.read('Catalogues/vandokkum15.csv',delimiter=',',format='no_header')
vcoords=SC(ra=data['col2'].data,dec=data['col3'].data)
plt.plot(vcoords.ra,vcoords.dec,'o',label='vanDokkum+15')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.legend(loc='best',ncol=2)
plt.show()

