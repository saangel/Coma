import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as p
import subprocess as s
import glob as g
import os

redshift = 0.0231 # redshift of the Coma Cluster
home=os.getcwd()
os.chdir("/home/simon/Documents/PEGASE-HR/bin")
prefix_run="test"
s.call("rm "+prefix_run+"*.dat",shell=True)
name_SSPs_input="SSPs_input"
s.call(["rm",name_SSPs_input])
SSPs_input=open(name_SSPs_input,"w")
SSPs_input.write("2\n0.08\n120\nA\ny\n1\n"+prefix_run) # IMF - lower mass - higher mass - SN model - Stellar wind - Spectral library
SSPs_input.close()
s.call("./SSPs_HR <"+name_SSPs_input,shell=True)
tracks=g.glob("*Z*.dat")
Zs=sorted([t[12:-4] for t in tracks])
scenario_name="scenario"
s.call(["rm",scenario_name])
name_scenario_input="scenario_input"
s.call(["rm",name_scenario_input])
scenario_input=open(name_scenario_input,"w")
scenario_input.write(scenario_name+"\n")
scenario_input.write(prefix_run+"_SSPs.dat\n")
scenario_input.write("0.05\n1\n")
props="n0y0nn0" # no infall - single burst SFH - evolution of metallicity - no substellar objects - no winds/nebular emission - no extinction
spectra_list=[]
for i,z in enumerate(Zs):
    s.call(["rm","spectra_"+str(i)+".fits"])
    scenario_input.write("spectra_"+str(i)+".fits\n")
    spectra_list.append("spectra_"+str(i)+".fits")
    scenario_input.write(str(z)+"\n")
    dump=[scenario_input.write(pr+"\n") for pr in props]
scenario_input.write("end")
scenario_input.close()
s.call("./scenarios_HR <"+name_scenario_input,shell=True)
name_spectra_input="spectra_input"
s.call(["rm",name_spectra_input])
spectra_input=open(name_spectra_input,"w")
spectra_input.write(scenario_name)
spectra_input.close()
s.call("./spectra_HR <"+name_spectra_input,shell=True)
for sp in spectra_list:
    hdu=p.open(sp,mode="update")
    wl=hdu[3].data["BFIT"]
    hdu[3].data["BFIT"]=[(1+redshift)*w for w in wl]
    hdu.flush()
    hdu.close()
name_color_input="color_input"
for i,sp in enumerate(spectra_list):
    s.call(["rm",name_color_input])
    s.call(["rm","color_"+str(i)+".dat"])
    color_input=open(name_color_input,"w")
    color_input.write(sp+"\n")
    color_input.write("color_"+str(i)+".dat")
    color_input.close()
    s.call("./colors_HR <"+name_color_input,shell=True)

os.chdir(home)
