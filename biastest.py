from astropy.io import fits as p
from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord as SC
from astropy import units as u
from astropy.wcs import WCS
from glob import glob as g
from sklearn.neighbors import KernelDensity as KD
from string import ascii_uppercase as LETTERS
import numpy as np
import os

'''
def makeslices(data):
    dataslices=[]
    for s in np.array_split(data,2): dataslices.append(np.array_split(s,2,axis=1))
    dataslices=[dataslices[0][0],dataslices[0][1],dataslices[1][0],dataslices[1][1]]
    return dataslices
def plotslices(dataslices,kdensity,rootshow,redrange,offsets,log):
    offset=[]
    for i,s in enumerate(dataslices):
        X=s.flatten()
        X=X[(np.isfinite(X)) & (X>0)]
        plt.xscale("log");plt.yscale("log")
        if kdensity:
            X=X[:,np.newaxis]
            x_plot=np.logspace(-9,np.log10(25),5000)[:,np.newaxis]
            kde=KD(kernel="epanechnikov",bandwidth=0.005).fit(X)
            print "fitted KDE"
            ld=kde.score_samples(x_plot)
            plt.plot(x_plot[:,0],np.exp(ld),'-',label=str(i))
        else:
            if redrange:
                if log: bins=np.logspace(np.log10(0.0012),np.log10(0.023),1000)
                else: bins = np.linspace(0.00012,0.0023,1000)
            else: bins=np.logspace(np.log10(min(X)),np.log10(max(X)),100)
            N,bins,patches=plt.hist(X,bins=bins,histtype="step",label=str(i),lw=3,alpha=.5)
            if rootshow:
                fit=np.poly1d(np.polyfit(bins[:-1],N,5))
                plt.plot(bins[:-1],fit(bins[:-1]),"-",label=str(i)+" fit",c="gray")
        if rootshow:
            roots=fit.deriv().r
            root=roots[np.isreal(roots)]
            root=np.real(root[root>0][0])
            plt.plot(root,fit(root),"ko")
            offset.append(root)
    plt.legend(loc="best")
    plt.show()
    if offsets: return offset
'''

bias=p.open("bias.fits")
sci_orig=p.open("h_v22_F814W_ivm_drz_cl_ver2.fits")
sci_new=p.open("test_sci.fits")
pixscale=sci_new[0].header["D001SCAL"]
fullbias=np.vstack((bias[1].data,bias[4].data))
data=sci_new[0].data
data_survey=sci_orig[1].data
#plotslices(makeslices(data_survey),False,False,True,False,False)

xdiff=np.abs(fullbias.shape[0]-data.shape[0])
ydiff=np.abs(fullbias.shape[1]-data.shape[1])

fullbias=np.pad(fullbias,((xdiff/2,xdiff/2+1),(ydiff/2,ydiff/2)),mode="constant",constant_values=0)
#offsets=plotslices(makeslices(data),False,False,True,True,False)
#biaslices=makeslices(fullbias)
#plotslices(biaslices,False,False,False,False)

#newbias=[bs-o for bs,o in zip(biaslices,offsets)]
#fullnewbias=np.vstack((np.hstack((newbias[0],newbias[1])),np.hstack((newbias[2],newbias[3])) ))

#plt.matshow(np.log10(data+(fullnewbias-fullbias)))
#plt.show()

regionfile=open("ACS_quadrants.reg")
regioninfo=regionfile.readlines()[3:]
wcs=WCS("test_sci.fits")
regions=[]

def makeRegions(hdu,layer):
    wcs=WCS(hdu[layer].header)
    image_empty=np.zeros_like(hdu[layer].data)

    image_empty[:]=np.nan
    regions=[]
    conds=[]
    plt.clf()

    f,ax=plt.subplots(1,figsize=(10,10))
    for r in regioninfo:
        r0=r[4:-2].split(",")
        ra,dec,sx,sy=r0[0],r0[1],float(r0[2][:-1])/pixscale/2,float(r0[3][:-1])/pixscale/2
        c=SC(ra,dec,unit=(u.hourangle,u.deg))
        out=wcs.wcs_world2pix(c.ra.value,c.dec.value,1)
        out.append(sx);out.append(sy)
        regions.append(out)

    regions=np.array(regions)
    regions=regions[[3,2,0,1]]
    slices=[]
    for i,r in enumerate(regions):
        data=hdu[layer].data
        sll=data[int(r[0]-r[2]):int(r[0]+r[2]),int(r[1]-r[3]):int(r[1]+r[3])]
        image_empty[int(r[0]-r[2]):int(r[0]+r[2]),int(r[1]-r[3]):int(r[1]+r[3])]=sll
        slices.append(sll)
        ax.annotate(str(i),xy=(r[0],r[1]),zorder=10)
    ax.matshow(data_survey,alpha=.7,cmap="gray",vmin=-0.025,vmax=0.25)
    ax.matshow(image_empty,vmin=-0.025,vmax=0.25)
#    plt.show()
    plt.savefig(hdu.fileinfo(layer)["filename"]+"_regions.png",density=300)
    plt.clf()
    return slices

dss=makeRegions(sci_orig,1)
ds=makeRegions(sci_new,0)
col=["black","red"]
typ=["survey","new"]

def plotting(j,ss):
    maxes=[]
    for i,a in enumerate(axarr.flat):
        s=ss[i]
        lim=.01
        s=s[s<lim]
        s=s[s>-lim]
        nbins=100
        bins=np.linspace(-lim,lim,nbins)
        xbins=np.linspace(-lim,lim,10000)
        n,bins,patches=a.hist(s.flatten(),bins=bins,histtype="step",label=typ[j],lw=3,alpha=.7,color=col[j])
        fit=np.poly1d(np.polyfit(bins[:-1],n,5))
        a.plot(bins[:-1],fit(bins[:-1]),"-",c="blue")
        argmax=np.argmax(fit(xbins))
        a.plot(xbins[argmax],fit(xbins[argmax]),"ko")
        maxes.append(xbins[argmax])
        a.set_title(str(i))
        a.axvline(0)
        a.set_xlim([-lim,lim])
        a.set_yscale("log")
        a.legend()
    return maxes
f, axarr = plt.subplots(2, 2,sharex=True,sharey=True,figsize=(10,8))
scimaxes=np.array(plotting(0,dss))
newmaxes=np.array(plotting(1,ds))
plt.savefig("slices_offsets",density=300)

diff=scimaxes-newmaxes
s=np.argmin(abs(diff))
diffs=diff-abs(diff[s])*(abs(diff[s])/diff[s])
print diffs
flcs=g("*_flc.fits")
'''
for f in flcs:
    j=0
    hdu=p.open(f)
    head=hdu[0].header
    exptime=head["EXPTIME"]
    fig,axarr=plt.subplots(2,3,figsize=(10,5))
    vmin=10
    vmax=40
    bins=np.linspace(vmin,vmax,100)
    xbins=np.linspace(vmin,vmax,10000)
    maxes=[]
    for i,a in enumerate(axarr[:,0]):
        im=a.matshow(hdu[1+3*i].data,vmin=vmin,vmax=vmax)
        splits=np.array_split(hdu[1+3*i].data,2,axis=1)
        for s in splits:
            n,bins,patches=axarr[1,1].hist(s.flatten(),label=str(j),bins=bins,histtype="step",lw=3,alpha=.5)
            fit=np.poly1d(np.polyfit(bins[:-1],n,9))
            axarr[1,1].plot(bins[:-1],fit(bins[:-1]),"-",c="black",alpha=.3)
            argmax=np.argmax(fit(xbins))
            axarr[1,1].scatter(xbins[argmax],fit(xbins[argmax]),marker="o",c="black")
            maxes.append(xbins[argmax])
            #leg=axarr[0,1].plot(j,xbins[argmax]-maxes[0],"ko",label="max")
            axarr[1,1].set_yscale("log")
            j=j+1
        axarr[1,1].legend()
    col=["red","green","blue"]
    plotlegend=[]
    props=["ATODGN","READNSE","BIASLEV"]
    pvals=np.empty((3,4))
    for i in range(4):
        for k,prop in enumerate(props):
            pvals[k,i]=head[prop+LETTERS[i]]
    axarr[0,1].scatter(maxes,pvals[0],label=props[0],s=15)
    axarr[0,2].scatter(maxes,pvals[1],label=props[1],s=15)
    axarr[1,2].scatter(maxes,pvals[2],label=props[2],s=15)

            #norm=1
            #if prop=="BIASLEV": norm=1000
            #t=axarr[0,1].scatter(i,head[prop+LETTERS[i]]-head[prop+LETTERS[0]],marker="o",color=col[k],label=prop)
            #plotlegend.append(t)
    axarr[0,1].legend(loc="best")
    axarr[0,1].set_xlabel("normalised quadrant maxes")
    axarr[0,1].set_ylabel("normalised chip property")
    #axarr[0,1].legend(handles=plotlegend[:3]+leg,loc="best")
    #axarr[0,1].set_yscale("symlog")
    plt.show()
'''
