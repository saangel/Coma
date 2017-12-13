import subprocess as s
import glob as g
import numpy as np
from astropy.io import fits as p
import os
home=os.getcwd()
workdir="/Users/simon/Documents/Coma/"
os.chdir(workdir)
cats=np.array(g.glob('Catalogues/*.dat'))

filters=['F336W','F475W','F814W','F160W']
table_names=[]
match_call="java -jar stilts.jar tmatchn multimode=group nin="+str(len(filters))+" matcher=sky params=0.5 "
for i,f in enumerate(filters):
	mask=np.array([f in c for c in cats])
	cat_subset=cats[mask]
	cat_subset=[c+"#2" for c in cat_subset]
	table_name="table_"+f+"_stilts.fits"
	print "joining "+table_name
	s.call("java -jar stilts.jar tcat in=' "+" ".join(cat_subset)+"' out="+table_name,shell=True)
	print "cleaning "+table_name
	table=p.open(table_name,mode="update")
	table_names.append(table_name)
	table[1].data=table[1].data[table[1].data["DELTA_J2000"]>27.9]
	table[1].data=table[1].data[table[1].data["DELTA_J2000"]<28.1]
	table.flush()
	table.close()
	match_call=match_call+"in"+str(i+1)+"="+table_name+" values"+str(i+1)+"='ALPHA_J2000 DELTA_J2000' "
match_call=match_call+"out=Coma_matched.fits"
print match_call
s.call(match_call,shell=True)
