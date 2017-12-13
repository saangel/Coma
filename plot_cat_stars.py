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
maxerr=.3
data_cat=data_cat[data_cat["MAGERR_APER_1"][:,1]<maxerr]
data_cat=data_cat[data_cat["MAGERR_APER_2"][:,1]<maxerr]
data_cat=data_cat[data_cat["MAGERR_APER_3"][:,1]<maxerr]
data_cat=data_cat[data_cat["MAGERR_APER_4"][:,1]<maxerr]
hdu.close()
hdu=p.open("../../Coding/synth/table_master.fits")
data_star=hdu[1].data
hdu.close()
ax1=plt.subplot2grid((8,10),(0,0),rowspan=4,colspan=4)
ax2=plt.subplot2grid((8,10),(0,5),sharex=ax1,sharey=ax1,rowspan=4,colspan=4)
ax3=plt.subplot2grid((8,10),(4,0),sharex=ax1,sharey=ax1,rowspan=4,colspan=4)
ax4=plt.subplot2grid((8,10),(4,5),sharex=ax1,sharey=ax1,rowspan=4,colspan=4)
ax1.set_xlim([-1,5])
ax1.set_ylim([-0.8,4])
axes=[ax1,ax2,ax3,ax4]
pos=[(0,4),(0,9),(4,4),(4,9)]
props=data_star.names[-4:]
vmins=(3000,-0.5,-4,-0.4)
vmaxs=(9000,6.5,1,1.21)
for pr,ax,ps,vmin,vmax in zip(props,axes,pos,vmins,vmaxs):
    ax.scatter(data_cat['MAG_APER_1'][:,1]-derred[0]-(data_cat["MAG_APER_3"][:,1]-derred[2]),data_cat["MAG_APER_3"][:,1]-derred[2]-(data_cat["MAG_APER_4"][:,1]-derred[3]),marker="o",c="black",edgecolors="none",s=5,zorder=2)
    fig=ax.scatter(data_star["f336w"]-data_star["f814w"],data_star["f814w"]-data_star["f160w"],c=data_star[pr],edgecolors="none",vmin=vmin,vmax=vmax)
    #ax.set_title(pr)
    cax=plt.subplot2grid((8,10),ps,rowspan=4)
    cb=plt.colorbar(fig,cax=cax)
    cax.set_title(pr,fontsize=10)
    #cb.set_label(pr,fontsize=5)
ax1.axes.get_xaxis().set_visible(False)
ax2.axes.get_xaxis().set_visible(False)
ax2.axes.get_yaxis().set_visible(False)
ax3.set_xlabel(r"F336W-F814W",fontsize=20)
ax3.set_ylabel(r"F814W-F160W",fontsize=20)
ax4.axes.get_yaxis().set_visible(False)
plt.show()
