#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from numpy import linalg as LA
from scipy.ndimage.filters import gaussian_filter
from util import taper_cosine, grad_struct, grad_cg

#Input paras
fdir        =sys.argv[1]
lambda0     =float(sys.argv[2])
eta         =float(sys.argv[3])
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
wdir        =optim['wdir']
fmod        =optim['fmod']
fmodref     =optim['fmodref']
fmod2       =optim['fmod2']
fgrad       =optim['fgrad']
fobstime    =optim['fobstime']
n2          =optim['nx']
n1          =optim['nz']
d2          =optim['dx']
d1          =optim['dz']
m2          =optim['mx']
m1          =optim['mz']
n2s         =optim['nxs']
n2e         =optim['nxe']
n1s         =optim['nzs']
n1e         =optim['nze']
m2s         =optim['mxs']
m2e         =optim['mxe']
m1s         =optim['mzs']
m1e         =optim['mze']
smooth_size =optim['smooth_size']
maxmode     =optim['maxmode']

# interpretation
modes = maxmode + 1
fdiff = fobstime[:-4] + '_difforig.csv'
sigma=smooth_size/np.sqrt(8)/d2
print('sigma(xs, zs):',sigma)

# read observed data
df = pd.read_csv(fdiff)
num, ncols = df.shape
header = df.columns.to_list()
if maxmode == 0:
    Cs =[i for i in header if 'M0_T' in i]
else:
    Cs = [i for i in header if '_T' in i]
diff = np.array(df[Cs], np.float32)
diff = np.where(diff > -900.0, diff, 0.0)
evxs = df['evx']
stxs = df['stx']
nper = len(Cs)

# read reference model
m0 = np.fromfile(fmod, np.float32)
mref = np.fromfile(fmodref, np.float32)
mrho = np.fromfile(fmod2, np.float32)

# read phase velocity
phvall = np.zeros([n2, nper])
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

#find station index
Icoords = np.zeros([num, 2], np.float32)
for idx in range(num):
    evx, stx = evxs[idx], stxs[idx]
    isrc, irec = (np.abs(coordx1d - evx)).argmin(), (np.abs(coordx1d - stx)).argmin()
    Icoords[idx, :] = min(isrc, irec), max(isrc, irec)

#sum up time diff. for each phase
np.seterr(divide='ignore', invalid='ignore')
gradfinal = np.zeros([n2,n1], np.float32)
gradsurf = np.zeros([n2, nper], np.float32)
taper_width=4.0 #smooth radius is 4 elements
for iph in range(nper):
    for i in range(num):
        isrc, irec = Icoords[i, :]
        taper1 = taper_cosine(1.0, n2, isrc, irec, taper_width, taper_width)
        gradsurf[:, iph] += taper1*diff[i, iph]
    tmp = np.where(phvall[:, iph]>0.0, -2.0*d2*gradsurf[:, iph]/phvall[:, iph]**2, 0.0)
    gradsurf[:, iph]=gaussian_filter(tmp, sigma, mode='nearest')#sigma
    #read kernels
    fname = 'ker_%s.nc'%Cs[iph]
    fin = os.path.join(wdir, fname)
    ncin = Dataset(fin)
    dcdm = ncin.variables['vs'][:,:].data
    ncin.close()
    gradfinal[:] += np.einsum('i,ij->ij',gradsurf[:, iph], dcdm)

# sum up all terms
## read reference model
m0 = np.fromfile(fmod, np.float32).reshape([n2, n1])
mref = np.fromfile(fmodref, np.float32).reshape([n2, n1])
mrho = np.fromfile(fmod2, np.float32).reshape([m2, m1])
c2, c1 = m0.shape
m0 = m0.flatten()
mref = mref.flatten()
mrho = mrho.flatten()
g2 = grad_struct(m0, mref, c1, c2, d1, d2)

m0 = m0.reshape([n2, n1])[n2s:n2e+1,n1s:n1e+1]
mref = mref.reshape([n2, n1])[n2s:n2e+1,n1s:n1e+1]
mrho = mrho.reshape([m2, m1])[m2s:m2e+1,m1s:m1e+1]
c2, c1 = m0.shape
m0 = m0.flatten()
mref = mref.flatten()
mrho = mrho.flatten()
g3 = grad_cg(mrho, m0, c1, c2, d1, d2)
#grad2 = np.zeros([n2,n1], np.float32)
grad3 = np.zeros([n2,n1], np.float32)
grad2 = g2.reshape([n2, n1])
grad3[n2s:n2e+1,n1s:n1e+1] = g3.reshape([c2, c1])
gradfinal = gradfinal + lambda0 * grad2 + eta * grad3

######### preprocessing #########
where_are_NaNs=np.isnan(gradfinal)
gradfinal[where_are_NaNs] = 0.0 #unknown error
gradfinal[:]=gaussian_filter(gradfinal, 4.0/np.sqrt(8), mode='nearest')#sigma
##################################
gradfinal.tofile(fgrad)

print("Finished...", __file__)
