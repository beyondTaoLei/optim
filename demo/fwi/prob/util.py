#!/usr/bin/env python3
"""
useful functions
Usage:
    [from util import fun] in another *.py file
https://stackoverflow.com/questions/20309456/call-a-function-from-another-file
"""
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from scipy.signal import medfilt2d
from scipy.interpolate import RectBivariateSpline

def taper_cosine(dt, nt, win_left, win_right,taper_width, verbose=False):
    """
    win_left: left window boundary [s]
    win_right: right window boundary [s]
    taper_width: width of the cosine taper applied to the seismogram [s]
    verbose: additional information as screen output if 'True'
    
    taper_width should be less than half of (win_right-win_left)
    copied from ses3d software (Andreas Fichtner)
    """
    #- compute the cosine taper and apply to seismogram
    t=np.zeros(nt,dtype=float)
    for i in range(nt):
        t[i]=dt*i #dt*(i+1)
    N_0=int(round(win_left/dt))
    N=int(round(win_right/dt)-N_0)
    
    taper=np.zeros(nt,dtype=float)
    n=(t>win_left) & (t<win_right)
    taper[n]=1.0

    n=t<=taper_width
    h=0.5*(1.0-np.cos(np.pi*t[n]/taper_width))
    taper[N_0:N_0+len(h)]=h
    h=0.5*(1.0+np.cos(np.pi*t[n]/taper_width))
    taper[N_0+N-len(h):N_0+N]=h
    if taper[N_0-1]==1.0:
        taper[N_0]=0.0
    if taper[N_0+N]==1.0:
        taper[N_0+N]=0.0
        
    if verbose==True:
      plt.plot(t,taper)
      plt.grid(True)
      plt.xlabel('time [s]')
      plt.title('cosine taper')
      plt.show()
      
    return taper

def extrapolation2d(data, size):
    """
    extrapolates 2D array on the four boundaries
    """
    row,col=data.shape
    row2=row+2*size
    col2=col+2*size
    data_new=np.zeros((row2,col2),dtype=data.dtype)
    
    #copy data into data_new
    for i in range(row):
        ii=i+size
        for j in range(col):
            jj=j+size
            data_new[ii][jj]=data[i][j]
    
        
    #left/right boundary
    for ii in range(size,row+size):
        for jj in range(0,size,1):
            data_new[ii][jj]=data_new[ii][size]
        for jj in range(col+size,col2,1):
            data_new[ii][jj]=data_new[ii][col+size-1]
            
    #top/bottom boundary incl. corners
    for jj in range(col2):
        for ii in range(0,size,1):
            data_new[ii][jj]=data_new[size][jj]
        for ii in range(row+size,row2,1):
            data_new[ii][jj]=data_new[row+size-1][jj]
            
    return data_new

def extract2d(data,r1,r2,c1,c2):
    """
    extracts the part of 2d array
    """
    row,col=data.shape
    if (0<=r1 and r1<=r2 and r2<= row) and \
     (0<=c1 and c1<=c2 and c2<= col):
        row2=r2-r1+1
        col2=c2-c1+1
        data_new=np.zeros((row2,col2),dtype=data.dtype)
        for i in range(row2):
            ii=i+r1
            for j in range(col2):
                jj=j+c1
                data_new[i][j]=data[ii][jj]
        return data_new
    else:
        return None

def medfilt2d_enhanced(arr,win):
    """
    because function scipy.signal.medfilt2d() perform badly
    on the boundary, we should extend the array on the outside
    of the boundary.
    """
    if win%2!=0:
        size=win//2
    else:
        raise ValueError("Please set an odd size for smooth window")
    row,col=arr.shape
    arr2=extrapolation2d(arr,size)
    arr3=medfilt2d(arr2,kernel_size=win)
    arr4=extract2d(arr3,size,row+size-1,size,col+size-1)
    return arr4

def fun(a,b):
    return a+b
