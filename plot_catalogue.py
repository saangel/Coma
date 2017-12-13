import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as p
from matplotlib.colors import LogNorm

import subprocess as s
import glob as g
import os
hdu=p.open("../Catalogues/Coma_matched.fits")
derred=[0.0429,0.0312,0.0143,0.0048]
data_cat=hdu[1].data
print len(data_cat)
hdu.close()
hdu=p.open("../../Coding/synth/table_master.fits")
data_stars=hdu[1].data
hdu.close()
max_err=.3
data_cat=data_cat[data_cat["MAGERR_APER_1"][:,1]<max_err]
print len(data_cat)
data_cat=data_cat[data_cat["MAGERR_APER_2"][:,1]<max_err]
print len(data_cat)
data_cat=data_cat[data_cat["MAGERR_APER_3"][:,1]<max_err]
print len(data_cat)
data_cat=data_cat[data_cat["MAGERR_APER_4"][:,1]<max_err]
print len(data_cat)
# rms_fluxradius=np.sqrt(np.sum([np.power(data_cat["FLUX_RADIUS_"+str(i+1)],2) for i in range(3)],axis=0)/3)
# deltatot=np.sqrt(np.sum([np.power(data_cat["MAGERR_APER_"+str(i+1)],2) for i in range(4)],axis=0)/4)
# smtot=np.sqrt(np.sum([np.power(data_cat["SPREAD_MODEL_"+str(i+1)],2) for i in range(4)],axis=0)/4)
# fluxmaxtot=np.sqrt(np.sum([np.power(data_cat["FLUX_MAX_"+str(i+1)],2) for i in range(4)],axis=0)/4)

rms_fluxradius=np.median([data_cat["FLUX_RADIUS_"+str(i+1)] for i in range(4)],axis=0)
deltatot=np.median([data_cat["MAGERR_APER_"+str(i+1)][:,1] for i in range(4)],axis=0)
smtot=np.median([data_cat["SPREAD_MODEL_"+str(i+1)] for i in range(4)],axis=0)
fluxmaxtot=np.median([data_cat["FLUX_MAX_"+str(i+1)] for i in range(4)],axis=0)

ax1=plt.subplot2grid((12,14),(0,0),rowspan=6,colspan=6)
ax2=plt.subplot2grid((12,14),(0,7),sharex=ax1,sharey=ax1,rowspan=6,colspan=6)
ax3=plt.subplot2grid((12,14),(6,0),sharex=ax1,sharey=ax1,rowspan=6,colspan=6)
ax4=plt.subplot2grid((12,14),(6,7),sharex=ax1,sharey=ax1,rowspan=6,colspan=6)
fig1=ax1.scatter(data_cat['MAG_APER_1'][:,1]-derred[0]-(data_cat["MAG_APER_3"][:,1]-derred[2]),data_cat["MAG_APER_3"][:,1]-derred[2]-(data_cat["MAG_APER_4"][:,1]-derred[3]),marker="o",c=deltatot,edgecolors="none",s=5)
ax1.set_xlim([-1,5])
ax1.set_ylim([-0.8,4])
cax=plt.subplot2grid((12,14),(0,6),rowspan=6)
cb=plt.colorbar(fig1,cax=cax)
cax.set_title(r"$\Delta$mag",fontsize=10)
cb.ax.tick_params(direction='in')
ax1.axes.get_xaxis().set_visible(False)
fig2=ax2.scatter(data_cat['MAG_APER_1'][:,1]-derred[0]-(data_cat["MAG_APER_3"][:,1]-derred[2]),data_cat["MAG_APER_3"][:,1]-derred[2]-(data_cat["MAG_APER_4"][:,1]-derred[3]),marker="o",c=smtot,edgecolors="none",norm=LogNorm(vmax=0.1,vmin=0.0001),s=5)
cax=plt.subplot2grid((12,14),(0,13),rowspan=6)
cb=plt.colorbar(fig2,cax=cax)
cax.set_title("SPREAD_MODEL",fontsize=10)
cb.ax.tick_params(direction='in')
ax2.axes.get_xaxis().set_visible(False)
ax2.axes.get_yaxis().set_visible(False)
fig3=ax3.scatter(data_cat['MAG_APER_1'][:,1]-derred[0]-(data_cat["MAG_APER_3"][:,1]-derred[2]),data_cat["MAG_APER_3"][:,1]-derred[2]-(data_cat["MAG_APER_4"][:,1]-derred[3]),marker="o",c=fluxmaxtot,edgecolors="none",norm=LogNorm(vmax=100),s=5)
cax=plt.subplot2grid((12,14),(6,6),rowspan=6)
cb=plt.colorbar(fig3,cax=cax)
cax.set_title(r"FLUX_MAX (counts/s)",fontsize=10)
cb.ax.tick_params(direction='in')
ax3.set_xlabel(r"F336W-F814W",fontsize=20)
ax3.set_ylabel(r"F814W-F160W",fontsize=20)
fig4=ax4.scatter(data_cat['MAG_APER_1'][:,1]-derred[0]-(data_cat["MAG_APER_3"][:,1]-derred[2]),data_cat["MAG_APER_3"][:,1]-derred[2]-(data_cat["MAG_APER_4"][:,1]-derred[3]),marker="o",c=rms_fluxradius,edgecolors="none",norm=LogNorm(vmax=50),s=5)
cax=plt.subplot2grid((12,14),(6,13),rowspan=6)
cb=plt.colorbar(fig4,cax=cax)
cax.set_title(r"FLUX_RADIUS (pix)",fontsize=10)
cb.ax.tick_params(direction='in')
ax4.axes.get_yaxis().set_visible(False)
plt.show()
