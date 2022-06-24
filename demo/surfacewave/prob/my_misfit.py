#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from numpy import linalg as LA

#Input paras
fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
wdir        =optim['wdir']
fobstime    =optim['fobstime']
fcost       =optim['fcost']
nx          =optim['nx']
maxmode     =optim['maxmode']
modes       =maxmode + 1
fsyntime    =fobstime[:-4] + '_syn.csv'
fdiff       =fobstime[:-4] + '_difforig.csv'

# read observed data
dfobs = pd.read_csv(fobstime)
num, ncols = dfobs.shape
header = dfobs.columns.to_list()
Cs = [i for i in header if '_T' in i]
phtimeobs = np.array(dfobs[Cs], np.float32)
nper = len(Cs)
evxs = dfobs['evx']
stxs = dfobs['stx']

# read phase velocity
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
d = {'kevnm':list(dfobs['kevnm'])}
df = pd.DataFrame(data = d)
df['kstnm'] = dfobs['kstnm']
df['evx']  = dfobs['evx']
df['evy'] = dfobs['evy']
df['stx']  = dfobs['stx']
df['sty'] = dfobs['sty']
df['dist'] = dfobs['dist']
for i in range(len(Cs)):
    df[Cs[i]] = phtime[:,i]
df.to_csv(fsyntime,float_format="%.4f",index=False)

#time difference
diff=np.where((phtime>0.0) &(phtimeobs>0.0), phtime-phtimeobs, -1000.0)
for i in range(len(Cs)):
    df[Cs[i]] =diff[:,i]
print('per.     total     diff')
for cc in Cs:
    print("%10s %6d %6d"%(cc, df.shape[0], sum(abs(df[cc])<900.0)))

print('per.     total     diff')
#save the difference
df.to_csv(fdiff,float_format="%.4f",index=False)

# save to fcost file
diff_filtered=diff[diff>-900.0]
misfit=0.5*LA.norm(diff_filtered)**2
fp = open(fcost,'w')
fp.write('%.6e\n'% misfit)
fp.close()

print("Finished...", __file__)