#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
from netCDF4 import Dataset
sys.path.append('/home/tao/seis_software/ATT_surfwaves/forward')
from dpcxx import GradientEval

def fnm_comm(mode, per):
    return 'M%1d_T%.2f'%(mode, per)

#input paras
fdir        =sys.argv[1]
mode        =int(sys.argv[2]) #0, 1
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
wavetype    =optim['wavetype']
wdir        =optim['wdir']
fmod        =optim['fmod']
nz          =optim['nz']
dz          =optim['dz']
lastele     =3
# vs perturbation (km/s): avoid the homogeneous media
if 'pertb' in optim.keys():
    pertb = optim['pertb']
else:
    pertb = 0.001

# read phase velocity
fnm = os.path.join(wdir, 'vmap_M%1d.nc'%mode)
ncin = Dataset(fnm)
nx = ncin.dimensions['I'].size
nper = ncin.dimensions['K'].size
coordx = ncin.variables['coordx'][:,:].data
period = ncin.variables['period'][:,:].data
phvel = ncin.variables['vel'][:,:].data
ncin.close()
periods = period[0, :]
coordx1d = coordx[:,0]
freqs = 1.0/periods

# read model
pertbs = (-1.0)**np.arange(nz)*pertb
vp2d = np.fromfile(fmod[:-3]+'.vp', dtype=np.float32).reshape([nx, nz])
vs2d = np.fromfile(fmod[:-3]+'.vs', dtype=np.float32).reshape([nx, nz])
rho2d = np.fromfile(fmod[:-3]+'.rho', dtype=np.float32).reshape([nx, nz])
vp2d /=1000.0
vs2d /=1000.0
rho2d /=1000.0
vs2d += np.tile(pertbs,[nx, 1])

# run kernel generation
phvel[:] /= 1000.0 #km/s
model = np.zeros([nz, 5])
model[:,0] = np.arange(nz) + 1.0
coordz1d = np.arange(nz)*dz
model[:,1] = coordz1d/1000.0 #'kilometers'
kernels = np.zeros([nper, nx, nz], np.float32)
for i in range(nx):
    model[:,2] = rho2d[i,:]
    model[:,3] = vs2d[i,:]
    model[:,4] = vp2d[i,:]
    cs = phvel[i, :]
    gradEval = GradientEval(model, wavetype)
    for f1, c1 in zip(freqs, cs):
        if c1 > 0.0:
            iper=(np.abs(periods - 1.0/f1)).argmin()
            kernels[iper, i, :] = gradEval.compute(f1, c1)
            # reset the last two elements
            kernels[iper, i, -lastele:] =kernels[iper, i, -lastele]

# save kernels
for iper in range(nper):
    ker=kernels[iper, :, :]
    where_are_NaNs=np.isnan(ker)
    ker[where_are_NaNs] = 0.0 #unknown error
    count=np.sum(where_are_NaNs)
    ############################
    # open a netCDF file to write
    fname='ker_'+fnm_comm(mode, periods[iper])+'.nc'
    fout=os.path.join(wdir, fname)
    ncfile = Dataset(fout, 'w', format='NETCDF3_CLASSIC')#'NETCDF4'
    #create dimensions
    ncfile.createDimension('I', nx)
    ncfile.createDimension('K', nz)
    ncfile.typevel='phv'
    #define variables
    kervs2d   = ncfile.createVariable('vs','f4',('I','K'))
    coordxout = ncfile.createVariable('coordx','f4',('I','K'))
    coordzout = ncfile.createVariable('coordz','f4',('I','K'))
    #set values
    kervs2d[:]   =ker
    coordxout[:] =np.tile(coordx1d, (nz,1)).T
    coordzout[:] =np.tile(coordz1d, (nx,1))
    #close ncfile 
    ncfile.close()
    print('NaN percentage: %s %.6f'%(fname, count/(nx*nz)))
    ############################

print("Finished...", __file__)