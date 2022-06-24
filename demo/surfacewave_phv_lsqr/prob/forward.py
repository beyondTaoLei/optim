#!/usr/bin/env python3
import os
import sys
import shutil
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from disba import PhaseDispersion
from scipy.ndimage.filters import gaussian_filter

#input paras
fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
fmod        =optim['fmod']
fperiods    =optim['fperiods']
algo        =optim['algo'] #'dunkin'
wavetype    =optim['wavetype']#'love' or 'rayleigh'
maxmode     =optim['maxmode']
wdir        =optim['wdir']
nx          =optim['n2g']
nz          =optim['n1g']
dx          =optim['d2g']
dz          =optim['d1g']
ox          =optim['o2g']
modes       =maxmode + 1

# vs perturbation (km/s): avoid the homogeneous media
if 'pertb' in optim.keys():
    pertb = optim['pertb']
else:
    pertb = 0.001

# make directory
if os.path.exists(wdir):
    shutil.rmtree(wdir, ignore_errors=True)
os.makedirs(wdir,exist_ok=True)

# interpretation
df = pd.read_csv(fperiods)
periods = np.array(df['T0'])
nper = len(periods)

# read model
pertbs = (-1.0)**np.arange(nz)*pertb
vp2d = np.fromfile(fmod[:-3]+'.vp', dtype=np.float32).reshape([nx, nz])
vs2d = np.fromfile(fmod[:-3]+'.vs', dtype=np.float32).reshape([nx, nz])
rho2d = np.fromfile(fmod[:-3]+'.rho', dtype=np.float32).reshape([nx, nz])
vp2d /=1000.0
vs2d /=1000.0
rho2d /=1000.0
vs2d += np.tile(pertbs,[nx, 1])
coordz1d = np.arange(nz)*dz
z = coordz1d/1000.0 #'kilometers'
tn = np.append(np.diff(z), [0, ])

# run FD
model=np.zeros([nz, 4])
model[:,0] = tn
phvel=np.zeros([modes, nx, nper], np.float32) #period from short to long
for i in range(nx):
    model[:,1] = vp2d[i,:]
    model[:,2] = vs2d[i,:]
    model[:,3] = rho2d[i,:]
    pd0 = PhaseDispersion(*model.T, algorithm=algo)
    for mode in range(modes):
        try:
            cp = pd0(periods, mode, wave=wavetype)
            pers, cs = cp.period, cp.velocity
            nper0 = len(pers)
        except:
            #smooth model
            print("smoothing at index: "+str(i))
            model[:,2] = gaussian_filter(model[:,2], 6, mode='nearest')
            pd0 = PhaseDispersion(*model.T, algorithm=algo)
            #retry
            try:
                cp = pd0(periods, mode, wave=wavetype)
                pers, cs = cp.period, cp.velocity
                nper0 = len(pers)
            except:
                nper0 = 0
                print("Warning: an exception occurred at index: "+str(i))
        #if nper0 > 0:
        phvel[mode, i, :nper0] = cs

phvel[:] *= 1000.0 #m/s
coordx1d = np.arange(nx)*dx+ox #coordinates
for mode in range(modes):
    # open a netCDF file to write
    fnm = os.path.join(wdir, 'vmap_M%1d.nc'%mode)
    ncfile = Dataset(fnm, 'w', format='NETCDF3_CLASSIC')#'NETCDF4'
    #create dimensions
    ncfile.createDimension('I', nx)
    ncfile.createDimension('K', nper)
    ncfile.typevel='phv'
    #define variables
    coordxout = ncfile.createVariable('coordx','f4',('I','K'))
    periodout = ncfile.createVariable('period','f4',('I','K'))
    phv2d = ncfile.createVariable('vel',    'f4',('I','K'))
    #set values
    coordxout[:] =np.tile(coordx1d,(nper,1)).T
    periodout[:] =np.tile(periods,(nx,1))
    phv2d[:]   = phvel[mode,:,:]
    #close ncfile
    ncfile.close()

print("Finished...", __file__)