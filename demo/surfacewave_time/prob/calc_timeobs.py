#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
from netCDF4 import Dataset

#Input paras
wdir        ="step1"
fperiods    ="list/periods.csv"
fstn_pairs  ="list/stn_pairs.csv"
fobstime    ="list/obs_time.csv"
nx          =1300
maxmode     =1
modes       =maxmode + 1

#read periods
dfper=pd.read_csv(fperiods)
periods =np.array(dfper['T0'])
nphase = len(periods)
Cs =['M%1d_T%.2f'%(mode, per) for mode in range(modes) for per in periods]
header=["kevnm", "kstnm", "evx", "evy", "stx", "sty" ,"dist", "Tmin", "Tmax"]+Cs
#kevnm,kstnm,comp,evla,evlo,stla,stlo,dist
stn_pairs = pd.read_csv(fstn_pairs)
num, ncols=stn_pairs.shape
evxs = stn_pairs['evx']
stxs = stn_pairs['stx']

# read phase velocity
nper = modes*nphase
phvall = np.zeros([nx, nper])
i0 = 0
for mode in range(modes):
    key='M%1d'%mode
    pers = [float(i[4:]) for i in header if key in i] #'M0_T0.70'
    #phase velocity
    fnm = os.path.join(wdir, 'vmap_M%1d.nc'%mode)
    ncin = Dataset(fnm)
    n2 = ncin.dimensions['I'].size
    nphase = ncin.dimensions['K'].size
    coordx = ncin.variables['coordx'][:,:].data
    period = ncin.variables['period'][:,:].data
    phv    = ncin.variables['vel'][:,:].data
    ncin.close()
    coordx1d = coordx[:,0]
    periods = period[0, :]
    idx=[]
    for per in pers:
        i, = np.where(np.isclose(periods, per))
        idx.append(i[0])
    phvall[:,i0:i0+len(idx)]=phv[:, idx]
    i0=i0+len(idx)

# phase velocity to duration in each cell
np.seterr(divide='ignore', invalid='ignore')
segment = np.diff(coordx1d, axis=0, append=2*coordx1d[-1]-coordx1d[-2])
deltaT = np.zeros_like(phvall)
for iph in range(nper):
    deltaT[:,iph] = np.where(phvall[:,iph]>0.0, segment/phvall[:,iph], -1.0)
    
#####################################
#generate the synthetic phase time
phtime = np.zeros([num, nper], np.float32)
for idx in range(num):
    evx, stx = evxs[idx], stxs[idx]
    isrc, irec = (np.abs(coordx1d - evx)).argmin(), (np.abs(coordx1d - stx)).argmin()
    isrc, irec = min(isrc, irec), max(isrc, irec)
    for iph in range(nper):
        if np.sum(deltaT[isrc:irec,iph]>0.0)==irec-isrc:
            phtime[idx,iph] = deltaT[isrc:irec,iph].sum()
        else:
            phtime[idx,iph] = -1000.0

#save phase time
d = {'kevnm':list(stn_pairs['kevnm'])}
df = pd.DataFrame(data=d)
df['kstnm'] = stn_pairs['kstnm']
df['evx'] = stn_pairs['evx']
df['evy'] = stn_pairs['evy']
df['stx'] = stn_pairs['stx']
df['sty'] = stn_pairs['sty']
df['dist'] = stn_pairs['dist']
df['Tmin'] = stn_pairs['Tmin']
df['Tmax'] = stn_pairs['Tmax']
for i in range(len(Cs)):
    df[Cs[i]] = phtime[:,i]
df.to_csv(fobstime,float_format="%.4f",index=False)
########################################

print("Finished...", __file__)