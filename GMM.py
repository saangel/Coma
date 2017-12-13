import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits as p
from sklearn.mixture import GaussianMixture as GMM
from sklearn.neighbors import KernelDensity as KDE
from itertools import combinations as comb

catalog_name="../Catalogues/Coma_matched_compact_test.fits"
#catalog_name="../Catalogues/Coma_matched_compact.fits"
c_hdu=p.open(catalog_name)
cdata=c_hdu[1].data

def gauss_function(x, amp, x0, sigma):
    return amp * np.exp(-(x - x0) ** 2. / (2. * sigma ** 2.))

plt.clf()
filters="ugiH"

combs=[cc for cc in comb(filters,2)]
aics=[]
bics=[]
ncompsmax=2
for nc in range(2,ncompsmax+1):
	gmm = GMM(n_components=nc, covariance_type="full", tol=0.001)
	ax1=plt.subplot2grid((3,3),(0,0))
	ax2=plt.subplot2grid((3,3),(0,1))
	ax3=plt.subplot2grid((3,3),(0,2))
	ax4=plt.subplot2grid((3,3),(1,1))
	ax5=plt.subplot2grid((3,3),(1,2))
	ax6=plt.subplot2grid((3,3),(2,2))
	axes=[ax1,ax2,ax3,ax4,ax5,ax6]
	for i,c in enumerate(combs):
		colour=cdata[c[0]]-cdata[c[1]]
		gmm=gmm.fit(np.expand_dims(colour,1))
		gmm_x=np.sort(colour)
		gmm_y=np.exp(gmm.score_samples(gmm_x.reshape(-1,1)))
		kde=KDE(bandwidth=(gmm_x[-1]-gmm_x[0])/30,kernel="epanechnikov")
		kde.fit(colour[:,np.newaxis])
		#axes[i].hist(colour, bins=np.linspace(gmm_x[0],gmm_x[-1],20), normed=True, alpha=0.5, color="blue")
		axes[i].fill(gmm_x,np.exp(kde.score_samples(gmm_x[:,np.newaxis])),ec="blue",fc="blue",alpha=.5)
		axes[i].plot(gmm_x, gmm_y, color="black", lw=2, label="GMM")
		axes[i].set_xlabel(c[0]+" - "+c[1],weight="heavy")
		for m, cc, w in zip(gmm.means_.ravel(), gmm.covariances_.ravel(), gmm.weights_.ravel()):
			gauss = gauss_function(x=gmm_x, amp=1, x0=m, sigma=np.sqrt(cc))
			gauss = gauss / np.trapz(gauss, gmm_x) * w
			axes[i].plot(gmm_x,gauss,color="grey")
		bic=gmm.bic(np.expand_dims(colour,1))
		aic=gmm.aic(np.expand_dims(colour,1))
		aics.append(aic); bics.append(bic)
		#axes[i].set_title("BIC,AIC = "+str(bic)+", "+str(aic),fontsize=10) 
		
	plt.show()
plt.clf()
aics=np.reshape(aics,(ncompsmax,6))
bics=np.reshape(bics,(ncompsmax,6))
fig=plt.figure()
for i in range(6):
	fig.add_subplot(2,3,i+1)
	plt.plot(range(1,ncompsmax+1),bics[:,i],"red",label="BIC")
	plt.plot(range(1,ncompsmax+1),bics[:,i],"ro")
	plt.plot(range(1,ncompsmax+1),aics[:,i],"blue",label="AIC")
	plt.plot(range(1,ncompsmax+1),aics[:,i],"bo")
	plt.title(str(combs[i][0]) + " - "+str(combs[i][1]))
	plt.yscale("symlog")
	plt.legend()
plt.show()
