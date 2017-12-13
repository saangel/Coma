from astropy.io import fits as p
from scipy import ndimage
from scipy import misc
import numpy as np
import glob as g
import subprocess as s
import os, platform, copy
import warnings
import time
from joblib import Parallel,delayed

home=os.getcwd()
workdir="/raid/saangel/Coma/Data"
sextractor="sex"
det_filt="/raid/saangel/local/share/sextractor/gauss_1.5_3x3.conv"

os.chdir(workdir)
params1='NUMBER\nALPHA_J2000\nDELTA_J2000\nXWIN_IMAGE\nYWIN_IMAGE\nFLUX_RADIUS\nKRON_RADIUS\nMAG_APER\nMAG_AUTO\nELONGATION\nTHETA_IMAGE\nA_IMAGE\nB_IMAGE\nFLAGS\nFLUX_APER\nFLUXERR_APER\nSNR_WIN\nVIGNET(10,10)'
defparam1=open('default1.param','w')
defparam1.write(params1)
defparam1.close()


files=sorted(g.glob("final_drz_sci_imgs/*sci.fits")+g.glob("Archive/Coma_Archive/HST/*/*sci.fits"))
weights=sorted(g.glob("final_drz_wht_imgs/*rms.wht.fits")+g.glob("Archive/Coma_Archive/HST/*/*rms.wht.fits"))

iters=3
med_multiplier=1.1
skips=1

def magic(sexfile,w):
    hdu=p.open(sexfile)
    bname=sexfile.split("/")[-1].split(".")[0]
    head=hdu[0].header
    data_orig=hdu[0].data
    data_orig=np.nan_to_num(data_orig)
    #med_size=np.int(border/60)
    med_size=29
    if med_size%2==0: med_size+=1
    vmin=-2
    vmax=2.5
    dataprime=copy.deepcopy(data_orig)
    for i in range(iters):
        #dataprime=copy.deepcopy(data_orig)
        if i==0: sex_file=sexfile
        else: sex_file=bname+'_op%.0f.fits'%i
        checkimage_name=bname+"_checkop%.0f.fits" %i
        cat_name=bname+"_op%.0f.cat" %i
        sex_file_data=p.open(sex_file)[0].data
        s.call([sextractor,sex_file+"[0]",'-PARAMETERS_NAME','default1.param','-FILTER_NAME',det_filt,\
        '-CATALOG_NAME',cat_name,'-CATALOG_TYPE','FITS_LDAC','-VERBOSE_TYPE','QUIET','-CHECKIMAGE_TYPE','SEGMENTATION'\
        ,'-CHECKIMAGE_NAME', checkimage_name,'-WEIGHT_TYPE','MAP_WEIGHT','-WEIGHT_IMAGE',w,\
        "-DETECT_THRESH","1","-ANALYSIS_THRESH","1"])
        cat_2=p.open(cat_name)
        cat_data=cat_2[2].data
        head2=cat_2[2].header
        a_image=cat_data["A_IMAGE"]
        b_image=cat_data["B_IMAGE"]
        reff=cat_data["FLUX_RADIUS"]
        not_galaxy_condition=((reff<np.mean(reff)+3*np.std(reff))&(cat_data["MAG_AUTO"]<99)&(a_image>0)&(b_image>0))
        not_gal_like=cat_data[not_galaxy_condition]
        not_gal_nums=not_gal_like["NUMBER"]    
        data_seg=p.open(checkimage_name)[0].data
        #fig, axs = plt.subplots(1,3)
        #fig.set_figwidth(24)
        #fig.set_figheight(8)
      #  axs[0,0].imshow(np.log10(sex_file_data),vmin=vmin,vmax=vmax)
      #  axs[0,0].set_
        t0=time.time()
        for n in not_gal_nums:
            mask = data_seg==n 
            dataprime=np.where(~mask,dataprime,0)
           # print '{0}'.format(float(n)/len(not_gal_nums))
        #if i == 0:
        #    compact_original=data_orig-dataprime
        #    p.writeto("compact_original.fits",not_gal_like,head2,clobber=True,output_verify='ignore')
        t1=time.time()
        print "Done masking, time ",t1-t0
        #axs[0].imshow(np.log10(dataprime),vmin=vmin,vmax=vmax)
        #axs[0].set_title("Point sources masked")
        t0=time.time()
        data_median=ndimage.median_filter(dataprime,med_size)
        t1=time.time()
        print "Done median-ing, time ",t1-t0
        #axs[1].imshow(np.log10(data_median),vmin=vmin,vmax=vmax)
        #axs[1].set_title("PS masked & medianed, medsize %.0f" %med_size)
        new_data=data_orig-data_median
        #new_data=new_data+np.abs(np.min(new_data))+0.001
        new_data=np.nan_to_num(new_data)
        #print np.abs(np.min(new_data))
        seg_image=bname+'_op%.0f.fits'%(i+1)
        p.writeto(seg_image,new_data,head,overwrite=True,output_verify='ignore')
        print "Next image saved"
       # im = axs[2].imshow(np.log10(new_data),vmin=vmin,vmax=vmax)
        #axs[2].set_title("original - median")
       # fig.colorbar(im,cax=axs[3])
      #  plt.show()
       # plt.clf()
        med_size=np.int(np.floor(med_size*med_multiplier))
        if med_size%2==0: med_size+=1
njobs=1
Parallel(n_jobs=njobs,verbose=10)(delayed(magic)(s,w) for w,s in zip(weights,files))


os.chdir(home)
