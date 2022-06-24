#!/usr/bin/env python3
"""
split CCs
Usage:
    python *
"""
import os, csv, math
import numpy as np
import pandas as pd

# the interstation distance is at least 1.5 times the wavelength 
# to approximately satisfy the far-field approximation of 
# surface-wave propagation (Yao et al., 2011).

#Input paras
coords=[]
for i in range(1,68):
    coords.append([10000.+i*8000.,500.])
Scoord=np.array(coords)
ntr=len(coords)
stn_list=["sz%03d"%(i+1) for i in range(ntr)]

mindist     =16.0        #minimum distance between two stations R1-R2 (km)
#for Tmax
vmax        =3.5 #the average phase veloctiy is ok 
nwavesmin   =1.5
#for Tmin
vmin        =3.0 #the minimum phase velocity
nwavesmax   =20
if vmin*nwavesmax<=vmax*nwavesmin:
    raise ValueError("change the pars: vmin*nwavesmax<vmax*nwavesmin")

#OUTPUT
fstn_info  ='list/stationinfo.csv'
fstn_pairs ='list/stn_pairs.csv'
fperiods ='list/periods.csv'

#coordinate 2
#station name
d={'station':stn_list}
df=pd.DataFrame(data=d)
df["x"] =Scoord[:,0]
df["y"] =Scoord[:,1]
df.to_csv(fstn_info,index=False)

#write station pairs.
row_list=[["kevnm", "kstnm", "evx", "evy", "stx", "sty" ,"dist", "Tmin", "Tmax"]]
for i in range(len(stn_list)):
    for j in range(i, len(stn_list)):#stn1 is on the left of stn2
        dist=math.sqrt((Scoord[i,0]-Scoord[j,0])**2+(Scoord[i,1]-Scoord[j,1])**2)/1000.0 #km
        if dist>=mindist:
            Tmin, Tmax=dist/vmin/nwavesmax, dist/vmax/nwavesmin
            str1=list([stn_list[i], stn_list[j], Scoord[i,0], Scoord[i,1],\
                    Scoord[j,0], Scoord[j,1], '%.2f'%dist, '%.2f'%Tmin, '%.2f'%Tmax])
            row_list.append(str1)

with open(fstn_pairs, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(row_list)

#periods
phase = ['ph_%d'%i for i in range(19)]
T0 = ['%.2f'%i for i in range(3, 40,2)]
d={'phase':phase}
df=pd.DataFrame(data=d)
df["T0"] =T0
df.to_csv(fperiods,index=False)


print("Finished...", __file__)