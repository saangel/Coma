import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as p
from astropy.table import Table as tab
import subprocess as s
import glob as g
import os
home=os.getcwd()
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
data_stars=hdu[1].data
hdu.close()
grids=sorted(g.glob("../PegaseTracks/color*dat"))
colormap = plt.cm.viridis
colors = [colormap(i) for i in np.linspace(0.1, 0.9,len(grids))]
markers=['o', 'v', '^', '*','d','H','D']
skips=1
for k,g in enumerate(grids):
    lines=open(g,"r")
    t=[]
    ug=[]
    gi=[]
    ih=[]
    ui=[]
    for i,l in enumerate(lines):
        if i == 22:
            met=float(l.split()[-1])
        if i >= 517:
            j=l.split()
            t.append(float(j[0]))
            ug.append(float(j[1]))
            gi.append(float(j[2]))
            ih.append(float(j[3]))
            ui.append(float(j[4]))
    t=np.array(t)
    ui=np.array(ui)
    gi=np.array(gi)
    ug=np.array(ug)
    ih=np.array(ih)
    pegase_table=tab([t,ui,gi,ug,ih],names=("time","u-i","g-i","u-g","i-h"))
    pegase_table.write("pegase_z="+str(met)+".fits",format="fits",overwrite=True)
    fig=plt.scatter(ui[::skips],ih[::skips],c=np.log10(t*1000000.)[::skips],marker=markers[k],s=100,label=str(met),edgecolors="none",zorder=5,alpha=1)
    plt.plot(ui,ih,'--',color='black',alpha=.7,zorder=3)
plt.scatter(data_cat['MAG_APER_1'][:,1]-derred[0]-(data_cat["MAG_APER_3"][:,1]-derred[2]),data_cat["MAG_APER_3"][:,1]-derred[2]-(data_cat["MAG_APER_4"][:,1]-derred[3]),marker='.',c='black',zorder=4,alpha=.3,edgecolors="none")
#plt.plot(data_stars['f336w']-data_stars["f814w"],data_stars["f814w"]-data_stars["f160w"],'.',color='purple',zorder=2,alpha=.8,marker="*",markeredgecolor="none")
plt.tight_layout(pad=0.1,h_pad=10, w_pad=None, rect=None)
plt.tick_params(axis='both',which='both',labelsize=15)
cb=plt.colorbar(fig)
cb.set_label('log age/yr',fontsize=15)
cb.ax.tick_params(labelsize=15) 
plt.xlabel("F336W-F814W",fontsize=15)
plt.ylabel("F814W-F160W",fontsize=15)
plt.xlim([-1,5])
plt.ylim([-0.8,4])
plt.legend(bbox_to_anchor=(1.2,-0.2), loc='lower right',ncol=7, fancybox=True, shadow=True,fontsize=15)
plt.show()




os.chdir(home)
