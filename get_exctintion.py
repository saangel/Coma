import astropy.coordinates as coord
import astropy.units as u
import glob as g
import numpy as np
from extinction import fitzpatrick99 as f99
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 
from astroquery.irsa_dust import IrsaDust
from astropy.io import fits as p



#coords=["02h25m59.00s    -04d29m40.00s","10h00m28.00s    02d12m30.00s","14h19m27.00s    52d40m56.00s","22h15m31.00s    -17d43m56.00s"]
coords=[]

hdu=p.open("../Catalogues/Coma_matched.fits")
data_cat=hdu[1].data
hdu.close()

ra_max=np.max(data_cat["ALPHA_J2000_1"][~np.isnan(data_cat["ALPHA_J2000_1"])])
ra_min=np.min(data_cat["ALPHA_J2000_1"][~np.isnan(data_cat["ALPHA_J2000_1"])])
dec_max=np.max(data_cat["DELTA_J2000_1"][~np.isnan(data_cat["DELTA_J2000_1"])])
dec_min=np.min(data_cat["DELTA_J2000_1"][~np.isnan(data_cat["DELTA_J2000_1"])])
ra_slices=20
dec_slices=20
for i,ii in enumerate(np.linspace(ra_min,ra_max,ra_slices)):
    for j,jj in enumerate(np.linspace(dec_min,dec_max,dec_slices)):
        coords.append([ii,jj])
        


'''filter_names="ugrizJHK"
filters=[]
for f in filter_names:
    filters.append((g.glob("synth/*"+f+"*Mega*")+g.glob("synth/"+f+"wircam*"))[0])
'''
filters=sorted(g.glob("../../Coding/synth/F*W.fil.txt"))
filtercolor=cm.jet(np.linspace(0,1,len(filters)))
#print filters
wl=np.empty(0)
filts=[]
wls=[]
ax1=plt.subplot2grid((3,1),(0,0))
ax2=plt.subplot2grid((3,1),(1,0),sharex=ax1)
ax3=plt.subplot2grid((3,1),(2,0))
ax1.set_xlim([2900,17100])
ax1.set_ylim([0,0.6])
spec=0.5
for i,f in enumerate(filters):
    wlr,eff=np.loadtxt(f,unpack=True)
    wls.append(wlr)
    wl=np.concatenate((wl,wlr))
    filts.append(eff)
    ax1.fill(wlr,eff,label=f.split('.')[0],edgecolor="none",color=filtercolor[i])
ax1.axhline(spec,color="black",lw=3,alpha=.5)
#    ax1.set_xlabel(r"$\lambda$ in $\AA$")
ax1.set_ylabel("Throughput")
ax1.axes.get_xaxis().set_visible(False)
wl=np.sort(wl)

corrections=np.empty((len(filters),len(coords)))
mags_notred=np.empty(len(filters))
mags_red=np.empty((len(filters),len(coords)))
alambdas=[ [[] for _ in coords] for _ in filts]
color=cm.viridis(np.linspace(0,1,len(coords)))
for i,c in enumerate(coords):
  C = coord.SkyCoord(str(c[0])+" "+str(c[1]),unit="deg",frame="fk5")
  table=IrsaDust.get_query_table(C,radius=None)
  eb_v=table["ext SandF mean"]
  #print eb_v.data[0]
  al_plot=f99(wl,eb_v.data[0]*3.1)
  for j,f in enumerate(filts):
      alambdas[j][i]=f99(wls[j],eb_v.data[0]*3.1)
  ax2.plot(wl,al_plot,label=str(c[0])[:6]+" "+str(c[1])[:4],color=color[i])
ax2.set_xlabel(r"$\lambda$ in $\rm \AA$")
ax2.set_ylabel("Extinction in magnitudes")
ax2.set_ylim([0,0.07])
alambdas=np.array(alambdas)

for j,f in enumerate(filts):
    diffs=np.gradient(wls[j])
    flux=sum(wls[j]*spec*f*diffs) #integration
    norm=sum(f*diffs/wls[j]) #normalisation following GALEXEV docs.
    for k,c in enumerate(coords):
        tau=alambdas[j,k]/2.5/np.log10(np.e)
        spec_red=spec*np.exp(-tau)
        flux_red=sum(wls[j]*spec_red*f*diffs)
        mags_red[j,k]=-2.5*np.log10(flux_red/norm)
    mags_notred[j]=-2.5*np.log10(flux/norm)


derred=[]
for j in range(len(coords)):
  #  print mags_notred-mags_red[:,j]
    derred.append(mags_notred-mags_red[:,j])
np.savetxt("derredening_indexes.txt",derred)
#ax2.legend(loc="best",ncol=ra_slices,fontsize=7)
#title=" ".join([filters[i].split(".")[0][6:]+" = "+str(m)[:8] for i,m in enumerate(mags_notred)])
#title=title+"\n"
#title=title+" ".join([filters[i].split(".")[0][6:]+" = "+str(m)[:8] for i,m in enumerate(mags_red[:,0])])
#plt.title(title)


derred=np.array(derred)

for i in range(4):
	ax3.hist(derred[:,i],label=filters[i][19:24],color=filtercolor[i])
ax3.legend(loc="best")
ax3.set_xlabel("delta mag")
ax3.set_ylabel("n")
plt.savefig("ExtinctionByBand.png",density=300)
plt.clf()
#for i,c in enumerate(coords):
	#for j in range(len(filters)):
		##print derred[i,j]
		#plt.scatter(c[0],c[1],s=50000*np.abs(derred[i,j]),facecolors="none",edgecolors=filtercolor[j])
#ax = plt.gca()
#ax.invert_xaxis()
#plt.xlabel("ra")
#plt.ylabel("dec")
#plt.show()


ax1=plt.subplot2grid((2,2),(0,0))
ax2=plt.subplot2grid((2,2),(0,1),sharex=ax1,sharey=ax1)
ax3=plt.subplot2grid((2,2),(1,0),sharex=ax1,sharey=ax1)
ax4=plt.subplot2grid((2,2),(1,1),sharex=ax1,sharey=ax1)

axes=[ax1,ax2,ax3,ax4]

for i,ax in enumerate(axes):
	ax.imshow(derred[:,i].reshape((ra_slices,dec_slices)),extent=[ra_max,ra_min,dec_min,dec_max])
	ax.plot(data_cat["ALPHA_J2000_4"],data_cat["DELTA_J2000_4"],"k,",alpha=1)
	ax.set_title(filters[i][19:24]+" "+str(np.min(derred[:,i]))[:8]+" - "+str(np.max(derred[:,i]))[:8])
plt.savefig("ExtinctionMap.png",density=300)


