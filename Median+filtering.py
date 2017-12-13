
# coding: utf-8

# In[9]:


from astropy.io import fits as p
from matplotlib import pyplot as plt
from scipy import ndimage
from scipy import misc
import numpy as np
import glob as g
import subprocess as s
import os, platform, copy
import warnings
import time
warnings.filterwarnings('ignore')
get_ipython().magic(u'matplotlib inline')


# In[10]:


home=os.getcwd()
workdir="/home/simon/Documents/Coma/First_Batch/"
sextractor="sextractor"
det_filt='/usr/share/sextractor/gauss_1.5_3x3.conv'
if platform.system() == "Darwin":
	workdir="/Users/simon/Documents/Coma/Data"
	sextractor="sex"
	det_filt="/usr/local/scisoft/packages/sextractor-2.5.0/config/gauss_1.5_3x3.conv"
os.chdir(workdir)
params1='NUMBER\nALPHA_J2000\nDELTA_J2000\nXWIN_IMAGE\nYWIN_IMAGE\nFLUX_RADIUS\nKRON_RADIUS\nMAG_APER\nMAG_AUTO\nELONGATION\nTHETA_IMAGE\nA_IMAGE\nB_IMAGE\nFLAGS\nFLUX_APER\nFLUXERR_APER\nSNR_WIN\nVIGNET(10,10)'
defparam1=open('default1.param','w')
defparam1.write(params1)
defparam1.close()
params2='NUMBER\nALPHA_J2000\nDELTA_J2000\nFLUX_RADIUS\nFWHM_IMAGE\nELLIPTICITY\nFWHMPSF_IMAGE\nMAG_AUTO\nMAGERR_AUTO\nMAG_PSF\nMAGERR_PSF\nMAG_MODEL\nMAGERR_MODEL\nSPREAD_MODEL\nSPREADERR_MODEL\nFLUX_MAX\nA_IMAGE\nB_IMAGE'
defparam2=open('default2.param','w')
defparam2.write(params2)
defparam2.close()


# In[11]:


files=sorted(g.glob("final_drz_sci_imgs/*160*sci.fits"))#+g.glob("*336*sci.fits")+g.glob("*475*sci.fits")+g.glob("*814*sci.fits")
weights=sorted(g.glob("final_drz_wht_imgs/*160*wht.fits"))#+g.glob("*336*sci.fits")+g.glob("*475*sci.fits")+g.glob("*814*sci.fits")


# In[12]:


#for f,w in zip(files,weights): print f, w


# In[13]:


conditions={"F336W": 1.8,"F475W": 1.9,"F814W": 1.9,"F160W": 2.15}


# In[14]:


conditions2={"F336W": 10,"F475W": 20,"F814W": 30,"F160W": 20}


# In[15]:


sexfile=files[0]
iters=3
med_multiplier=1.1
skips=40


# In[16]:


for j,sexfile in enumerate(files[::skips]):
    conforme=False
    hdu=p.open(sexfile)
    head=hdu[0].header
    data_orig=hdu[0].data
    data_orig=np.nan_to_num(data_orig)
    pix_scale=head["D001SCAL"]
    im_size=head['NAXIS1'],head['NAXIS2']
    zp=-2.5*np.log10(float(head['PHOTFLAM']))+float(head['PHOTZPT'])-5*np.log10(float(head['PHOTPLAM']))+18.6921
    i=0
    filt=head['PHOTMODE'].split()[2]
    border=min(im_size)
    #med_size=np.int(border/60)
    med_size=np.int(13*conditions[filt])
    if med_size%2==0: med_size+=1
    vmin=-2
    vmax=2.5
    w=weights[j]
    dataprime=copy.deepcopy(data_orig)
    s.call([sextractor,sexfile,'-PARAMETERS_NAME','default1.param','-FILTER_NAME',det_filt,                '-CATALOG_NAME',"original.cat",'-CATALOG_TYPE','FITS_LDAC','-VERBOSE_TYPE','QUIET',                '-WEIGHT_TYPE','MAP_WEIGHT','-WEIGHT_IMAGE',w,"-DETECT_THRESH","1.5","-ANALYSIS_THRESH","1.5"])
    for i in range(iters):
        #dataprime=copy.deepcopy(data_orig)
        if i==0: sex_file=sexfile
        else: sex_file='op%.0f.fits'%i
        checkimage_name="checkop%.0f.fits" %i
        cat_name="op%.0f.cat" %i
        sex_file_data=p.open(sex_file)[0].data
        t0=time.time()
        s.call([sextractor,sex_file,'-PARAMETERS_NAME','default1.param','-FILTER_NAME',det_filt,                '-CATALOG_NAME',cat_name,'-CATALOG_TYPE','FITS_LDAC','-VERBOSE_TYPE','QUIET',                '-CHECKIMAGE_TYPE','SEGMENTATION','-CHECKIMAGE_NAME', checkimage_name,                '-WEIGHT_TYPE','MAP_WEIGHT','-WEIGHT_IMAGE',w,"-DETECT_THRESH","1.5","-ANALYSIS_THRESH","1.5"])
        t1=time.time()
        print "Done SExtractor, time ",t1-t0
        cat_2=p.open(cat_name)
        cat_data=cat_2[2].data
        head2=cat_2[2].header
        a_image=cat_data["A_IMAGE"]
        b_image=cat_data["B_IMAGE"]
        not_galaxy_condition=((a_image<=3*conditions[filt])&(cat_data["MAG_AUTO"]<99)&(a_image>0)&(b_image>0))
        not_gal_like=cat_data[not_galaxy_condition]
        not_gal_nums=not_gal_like["NUMBER"]    
        data_seg=p.open(checkimage_name)[0].data
        fig, axs = plt.subplots(1,3)
        fig.set_figwidth(24)
        fig.set_figheight(8)
      #  axs[0,0].imshow(np.log10(sex_file_data),vmin=vmin,vmax=vmax)
      #  axs[0,0].set_title("Original image")
        t0=time.time()
        for n in not_gal_nums:
            mask = data_seg==n 
            dataprime=np.where(~mask,dataprime,0)
           # print '{0}'.format(float(n)/len(not_gal_nums))
        if i == 0:
            compact_original=data_orig-dataprime
            p.writeto("compact_original.fits",not_gal_like,head2,clobber=True,output_verify='ignore')
        t1=time.time()
        print "Done masking, time ",t1-t0
        axs[0].imshow(np.log10(dataprime),vmin=vmin,vmax=vmax)
        axs[0].set_title("Point sources masked")
        t0=time.time()
        data_median=ndimage.median_filter(dataprime,med_size)
        t1=time.time()
        print "Done median-ing, time ",t1-t0
        axs[1].imshow(np.log10(data_median),vmin=vmin,vmax=vmax)
        axs[1].set_title("PS masked & medianed, medsize %.0f" %med_size)
        new_data=data_orig-data_median
        #new_data=new_data+np.abs(np.min(new_data))+0.001
        new_data=np.nan_to_num(new_data)
        #print np.abs(np.min(new_data))
        seg_image='op%.0f.fits'%(i+1)
        p.writeto(seg_image,new_data,head,clobber=True,output_verify='ignore')
        print "Next image saved"
        im = axs[2].imshow(np.log10(new_data),vmin=vmin,vmax=vmax)
        axs[2].set_title("original - median")
       # fig.colorbar(im,cax=axs[3])
        plt.show()
        plt.clf()
        med_size=np.int(np.floor(med_size*med_multiplier))
        if med_size%2==0: med_size+=1
        #conf=raw_input(u"Iteración %.0f. Está conforme? (Y/N)\n" %(i+1))
        #if conf=="Y" or conf=="y": conforme=True
    p.writeto("compact_final.fits",data_orig-dataprime,head,clobber=True,output_verify='ignore')
    s.call([sextractor,sex_file,'-PARAMETERS_NAME','default1.param','-FILTER_NAME',det_filt,                '-CATALOG_NAME',"compact_final.cat",'-CATALOG_TYPE','FITS_LDAC','-VERBOSE_TYPE','QUIET',                '-WEIGHT_TYPE','MAP_WEIGHT','-WEIGHT_IMAGE',w,"-DETECT_THRESH","1.5","-ANALYSIS_THRESH","1.5"])
    params3='NUMBER\nALPHA_J2000\nDELTA_J2000\nFLUX_RADIUS\nELLIPTICITY\nMAG_APER(2)\nMAGERR_APER(2)\nFLUX_MAX\nA_IMAGE\nB_IMAGE'
    defparam3=open('default_comparison.param','w')
    defparam3.write(params3)
    defparam3.close()


# In[ ]:


h=p.open("compact_final.cat")
compact_data=h[2].data


# In[ ]:


compact_data.columns
cdata=compact_data[(compact_data["flux_radius"]<=4)&(compact_data["flux_radius"]>0)]


# In[ ]:


fig, axs = plt.subplots(1,1)
fig.set_figwidth(8)
fig.set_figheight(8)
plt.scatter(cdata["ALPHA_J2000"],cdata["DELTA_J2000"],c=cdata["flux_radius"])
plt.show()
plt.hist(cdata["flux_radius"],bins=30)
plt.show()
plt.plot(cdata["flux_radius"],cdata["MAG_Aper"],"o",alpha=.3)


# In[ ]:





# In[ ]:


s.call([sextractor,'compact_final.fits','-PARAMETERS_NAME','default_comparison.param','-FILTER_NAME',det_filt,                '-CATALOG_NAME',"compact_comparison.fits",'-CATALOG_TYPE','FITS_LDAC','-VERBOSE_TYPE','QUIET',                '-WEIGHT_TYPE','MAP_WEIGHT','-WEIGHT_IMAGE',w,"-DETECT_THRESH","1","-ANALYSIS_THRESH","1",           "-PHOT_APERTURES","4,6"])
s.call([sextractor,'compact_final.fits',sexfile,'-PARAMETERS_NAME','default_comparison.param','-FILTER_NAME',det_filt,                '-CATALOG_NAME',"compact_comparison_doubleimage.fits",'-CATALOG_TYPE','FITS_LDAC','-VERBOSE_TYPE','QUIET',                '-WEIGHT_TYPE','MAP_WEIGHT,MAP_WEIGHT','-WEIGHT_IMAGE',w+","+w,"-DETECT_THRESH","1","-ANALYSIS_THRESH","1",           "-PHOT_APERTURES","4,6"])


# In[ ]:


fig=plt.figure(figsize=(40,40))
plt.imshow(np.log10(data_orig),cmap="gray",alpha=1)
plt.imshow(np.log10(data_orig-dataprime),cmap="inferno",alpha=1)
#plt.show()
plt.imshow(np.log10(compact_original),cmap="viridis",alpha=1)
plt.show()


# In[37]:


vmin=-2
vmax=2.5
hdu=p.open(sexfile)
head=hdu[0].header
img=np.nan_to_num(hdu[0].data)
f = np.fft.fft2(img)
fshift = np.fft.fftshift(f)
magnitude_spectrum = np.log(np.abs(fshift))
fig=plt.figure(figsize=(20,20))
plt.subplot(121),plt.imshow(np.log10(img),vmin=vmin,vmax=vmax)
plt.title('Input Image'), plt.xticks([]), plt.yticks([])
plt.subplot(122),plt.imshow(magnitude_spectrum)
plt.title('Magnitude Spectrum'), plt.xticks([]), plt.yticks([])
plt.show()
rows, cols = img.shape
crow,ccol = rows/2 , cols/2
mask_size=1000

fshiftprime=copy.deepcopy(fshift)
fshiftprime[crow-mask_size:crow+mask_size, ccol-mask_size:ccol+mask_size] = 0
fshift = fshift - fshiftprime

#fshift[crow-mask_size:crow+mask_size, ccol-mask_size:ccol+mask_size] = 0

f_ishift = np.fft.ifftshift(fshift)
img_back = np.fft.ifft2(f_ishift)
img_back = np.abs(img_back)
fig=plt.figure(figsize=(20,20))
plt.subplot(131),plt.imshow(np.log10(img),vmin=vmin,vmax=vmax)
plt.title('Input Image'), plt.xticks([]), plt.yticks([])
plt.subplot(132),plt.imshow(np.log(np.abs(f_ishift)))
plt.title('HPF'), plt.xticks([]), plt.yticks([])
plt.subplot(133),plt.imshow(np.log10(img_back),vmin=vmin,vmax=vmax)
plt.title('Result in VIR'), plt.xticks([]), plt.yticks([])

plt.show()


# In[17]:


print img_back


# In[ ]:





# In[ ]:


print np.unique(new_data), len(np.unique(new_data))
plt.hist(new_data.flatten(),log=True,bins=100)
plt.show()
plt.hist(np.log10(new_data.flatten()),log=True,bins=100)
plt.show()
plt.hist(np.log10(new_data).flatten(),log=True,bins=100)
plt.show()


# In[ ]:


data_orig2=data_orig+np.abs(np.min(data_orig))+0.001
data_median2=data_median+np.abs(np.min(data_median))+0.001
plt.imshow(data_orig2)
plt.show()
plt.imshow(data_median2)
plt.show()
plt.hist(data_orig2.flatten(),log=True,bins=20)
plt.show()
plt.hist(data_median2.flatten(),log=True,bins=20)
plt.show()
plt.hist(np.log10(data_orig2.flatten()),log=True,bins=20,normed=True,histtype="step")
plt.hist(np.log10(data_median2.flatten()),log=True,bins=20,normed=True,histtype="step")
plt.show()


# In[ ]:


data_fft=np.fft.fft2(data_orig)


# In[ ]:


plt.figure(figsize=(20,20))
plt.imshow(np.log(np.abs(data_fft)))
plt.colorbar()
plt.show()
plt.figure(figsize=(20,20))
plt.imshow(np.log10(np.real(data_fft)))
plt.colorbar()
plt.show()
plt.figure(figsize=(20,20))
plt.imshow(np.log10(np.imag(data_fft)))
plt.colorbar()
plt.show()


# In[ ]:


print pix_scale


# In[ ]:


print head["PHOTZPT"]


# In[ ]:


print mask.shape
print dataprime.shape


# In[ ]:


t0= time.time()
for i in range(100000000):
    pass
print time.time()-t0


# In[ ]:


a=np.floor(2.222)
print a


# In[ ]:


conforme=False
hdu=p.open(sexfile)
head=hdu[0].header
data_orig=hdu[0].data
data_orig=np.nan_to_num(data_orig)
pix_scale=head["D001SCAL"]
im_size=head['NAXIS1'],head['NAXIS2']
zp=-2.5*np.log10(float(head['PHOTFLAM']))+float(head['PHOTZPT'])-5*np.log10(float(head['PHOTPLAM']))+18.6921
i=0
border=min(im_size)
a_limit=0.1
if "160" in sexfile: a_limit=0.8
med_size=np.int(3*a_limit/pix_scale)
if med_size%2==0: med_size=med_size+1
vmin=-2
vmax=2.5
iters=5
for i in range(iters):
    dataprime=copy.deepcopy(data_orig)
    if i==0: sex_file=sexfile
    else: sex_file='op%.0f.fits'%i
    checkimage_name="checkop%.0f.fits" %i
    cat_name_psf="op_psf%.0f.cat" %i
    psf_name="op_psf%.0f.psf" %i
    cat_name="op%.0f.cat" %i
    sex_file_data=p.open(sex_file)[0].data
    t0=time.time()
    s.call([sextractor,sex_file,'-PARAMETERS_NAME','default1.param','-FILTER_NAME',filt,            '-CATALOG_NAME',cat_name_psf,'-CATALOG_TYPE','FITS_LDAC','-VERBOSE_TYPE','QUIET',            "-DETECT_THRESH","1.5","-ANALYSIS_THRESH","1.5"])
    t1=time.time()
    print "Done SExtractor 1, time ",t1-t0
    t0=time.time()
    s.call(['psfex',cat_name_psf,'-PHOTFLUX_KEY','FLUX_APER','-PHOTFLUXERR_KEY','FLUXERR_APER'            ,'-PSFVAR_NSNAP','4','-PSF_SIZE','15,15','-CHECKPLOT_TYPE','NONE','-CHECKIMAGE_TYPE','NONE'            ,'-NTHREADS','0',"-CENTER_KEYS","XWIN_IMAGE,YWIN_IMAGE","-PSFVAR_KEYS","XWIN_IMAGE,YWIN_IMAGE"])
    t1=time.time()
    print "Done PSFEx, time ",t1-t0

    t0=time.time()
    s.call([sextractor,sex_file,'-PARAMETERS_NAME','default2.param','-FILTER_NAME',filt,            '-CATALOG_NAME',cat_name,'-CATALOG_TYPE','FITS_LDAC','-VERBOSE_TYPE','QUIET',            '-CHECKIMAGE_TYPE','SEGMENTATION','-CHECKIMAGE_NAME', checkimage_name,            "-DETECT_THRESH","1.5","-ANALYSIS_THRESH","1.5","-PSF_NAME",psf_name])
    t1=time.time()
    print "Done SExtractor 2, time ",t1-t0
    cat_data=p.open(cat_name)[2].data
    a_image=cat_data["A_IMAGE"]
    b_image=cat_data["B_IMAGE"]
   # not_galaxy_condition=(cat_data["SPREAD_MODEL"]<0.1) & (cat_data["SPREAD_MODEL"]>0)\
   # & (cat_data["FLUX_RADIUS"]>0) & (cat_data["FLUX_RADIUS"]<1.5*np.mean(cat_data["FLUX_RADIUS"]))
    not_galaxy_condition=(cat_data["SPREAD_MODEL"]<0.1) & (cat_data["SPREAD_MODEL"]>0)    & (cat_data["FLUX_RADIUS"]>0) & (cat_data["FLUX_MAX"]<5)  
    not_gal_like=cat_data[not_galaxy_condition]
    not_gal_nums=not_gal_like["NUMBER"]    
    print len(not_gal_like)
    data_seg=p.open(checkimage_name)[0].data
    fig, axs = plt.subplots(1,3)
    fig.set_figwidth(24)
    fig.set_figheight(8)
  #  axs[0,0].imshow(np.log10(sex_file_data),vmin=vmin,vmax=vmax)
  #  axs[0,0].set_title("Original image")
    t0=time.time()
    for n in not_gal_nums:
        mask = data_seg==n 
        dataprime=np.where(~mask,dataprime,0)
       # print '{0}'.format(float(n)/len(not_gal_nums))
    t1=time.time()
    print "Done masking, time ",t1-t0
    axs[0].imshow(np.log10(dataprime),vmin=vmin,vmax=vmax)
    axs[0].set_title("Point sources masked")
    t0=time.time()
    data_median=ndimage.median_filter(dataprime,med_size)
    t1=time.time()
    print "Done median-ing, time ",t1-t0
    axs[1].imshow(np.log10(data_median),vmin=vmin,vmax=vmax)
    axs[1].set_title("PS masked & medianed, medsize %.0f" %med_size)
    new_data=data_orig-data_median
    #new_data=new_data+np.abs(np.min(new_data))+0.001
    new_data=np.nan_to_num(new_data)
    #print np.abs(np.min(new_data))
    seg_image='op%.0f.fits'%(i+1)
    p.writeto(seg_image,new_data,head,clobber=True,output_verify='ignore')
    print "Next image saved"
    im = axs[2].imshow(np.log10(new_data),vmin=vmin,vmax=vmax)
    axs[2].set_title("original - median")
   # fig.colorbar(im,cax=axs[3])
    plt.show()
    plt.clf()
    med_size=np.int(np.floor(med_size*med_multiplier))
    if med_size%2==0: med_size=med_size+1
    #conf=raw_input(u"Iteración %.0f. Está conforme? (Y/N)\n" %(i+1))
    #if conf=="Y" or conf=="y": conforme=True
    #i=i+1


# In[ ]:


gv=(cat_data["SPREAD_MODEL"]>0)&(cat_data["FLUX_RADIUS"]>0)
vg=cat_data[gv]
plt.hist(vg["SPREAD_MODEL"],bins=200)
plt.xscale("log")
plt.show()


# In[ ]:


plt.plot(vg["SPREAD_MODEL"],vg["FLUX_RADIUS"],"o")
plt.xscale("log")
plt.yscale("log")
plt.show()


# In[ ]:




