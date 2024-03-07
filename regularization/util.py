#!/usr/bin/env python3
"""
useful functions
Usage:
    [from util import fun] in another *.py file
https://stackoverflow.com/questions/20309456/call-a-function-from-another-file
"""
import numpy as np
from scipy import sparse
from numpy import linalg as LA

def opt_identity(n1, n2, d1=1.0, d2=1.0):
    """
    for smallest model
    """
    n12 = n1*n2
    coeff = np.ones(n12, np.float32)
    W = sparse.diags(coeff, shape=(n12, n12))
    return W.tocsc()

def opt_deriv(n1, n2, d1=1.0, d2=1.0):
    """
    for flattest model
    """
    n12 = n1*n2
    coeff = np.ones(n12, np.float32)
    #n1 direction
    W1 = sparse.diags(-coeff, shape=(n12, n12))
    W1 += sparse.diags(coeff, offsets=1, shape=(n12, n12))
    W1 = W1.tolil()
    i = n1-1
    for j in range(n2):
        idx = j * n1 + i
        W1[idx,:] = W1[idx-1,:]
    
    #n2 direction
    W2 = sparse.diags(-coeff, shape=(n12, n12))
    W2 += sparse.diags(coeff, offsets=n1, shape=(n12, n12))
    W2 = W2.tolil()
    j = n2-1
    for i in range(n1):
        idx = j * n1 + i
        W2[idx,:] = W2[idx-n1,:]
    
    #different directions
    return (W1/d1).tocsc(), (W2/d2).tocsc()

def opt_laplace(n1, n2, d1=1.0, d2=1.0, a1=1.0, a2=1.0):
    """
    for smoothest model
    """
    n12 = n1*n2
    coeff = np.ones(n12, np.float32)
    #n1 direction
    W1 = sparse.diags(-2.0*coeff, shape=(n12, n12))
    W1 += sparse.diags(coeff, offsets=1, shape=(n12, n12))
    W1 += sparse.diags(coeff, offsets=-1, shape=(n12, n12))
    W1 = W1.tolil()
    i = 0
    for j in range(n2):
        idx = j * n1 + i
        W1[idx,:] = W1[idx+1,:]
    
    i = n1-1
    for j in range(n2):
        idx = j * n1 + i
        W1[idx,:] = W1[idx-1,:]
    
    #n2 direction
    W2 = sparse.diags(-2.0*coeff, shape=(n12, n12))
    W2 += sparse.diags(coeff, offsets=n1, shape=(n12, n12))
    W2 += sparse.diags(coeff, offsets=-n1, shape=(n12, n12))
    W2 = W2.tolil()
    j = 0
    for i in range(n1):
        idx = j * n1 + i
        W2[idx,:] = W2[idx+n1,:]
    
    j = n2-1
    for i in range(n1):
        idx = j * n1 + i
        W2[idx,:] = W2[idx-n1,:]
    
    #sum up
    W = a1 / d1 / d1 * W1 + a2 /d2 / d2 * W2
    
    return W.tocsc()

def opt_identity3d(n1, n2, n3, d1=1.0, d2=1.0, d3=1.0):
    """
    for smallest model
    """
    n123 = n1*n2*n3
    coeff = np.ones(n123, np.float32)
    W = sparse.diags(coeff, shape=(n123, n123))
    return W.tocsc()

def opt_deriv3d(n1, n2, n3, d1=1.0, d2=1.0, d3=1.0):
    """
    for flattest model
    """
    n123 = n1*n2*n3
    n12 = n1*n2
    coeff = np.ones(n123, np.float32)
    #n1 direction
    W1 = sparse.diags(-coeff, shape=(n123, n123))
    W1 += sparse.diags(coeff, offsets=1, shape=(n123, n123))
    W1 = W1.tolil()
    i = n1-1
    for k in range(n3):
        for j in range(n2):
            idx = k * n12 + j * n1 + i
            W1[idx,:] = W1[idx-1,:]
    
    #n2 direction
    W2 = sparse.diags(-coeff, shape=(n123, n123))
    W2 += sparse.diags(coeff, offsets=n1, shape=(n123, n123))
    W2 = W2.tolil()
    j = n2-1
    for k in range(n3):
        for i in range(n1):
            idx = k * n12 + j * n1 + i
            W2[idx,:] = W2[idx-n1,:]
    
    #n3 direction
    W3 = sparse.diags(-coeff, shape=(n123, n123))
    W3 += sparse.diags(coeff, offsets=n12, shape=(n123, n123))
    W3 = W3.tolil()
    k = n3-1
    for j in range(n2):
        for i in range(n1):
            idx = k * n12 + j * n1 + i
            W3[idx,:] = W3[idx-n12,:]
    
    #different directions
    return (W1/d1).tocsc(), (W2/d2).tocsc(), (W3/d3).tocsc()

def opt_laplace3d(n1, n2, n3, d1=1.0, d2=1.0, d3=1.0, a1=1.0, a2=1.0, a3=1.0):
    """
    for smoothest model
    """
    n123 = n1*n2*n3
    n12 = n1*n2
    coeff = np.ones(n123, np.float32)
    #n1 direction
    W1 = sparse.diags(-2.0*coeff, shape=(n123, n123))
    W1 += sparse.diags(coeff, offsets=1, shape=(n123, n123))
    W1 += sparse.diags(coeff, offsets=-1, shape=(n123, n123))
    W1 = W1.tolil()
    i = 0
    for k in range(n3):
        for j in range(n2):
            idx = k * n12 + j * n1 + i
            W1[idx,:] = W1[idx+1,:]
    
    i = n1-1
    for k in range(n3):
        for j in range(n2):
            idx = k * n12 + j * n1 + i
            W1[idx,:] = W1[idx-1,:]
    
    #n2 direction
    W2 = sparse.diags(-2.0*coeff, shape=(n123, n123))
    W2 += sparse.diags(coeff, offsets=n1, shape=(n123, n123))
    W2 += sparse.diags(coeff, offsets=-n1, shape=(n123, n123))
    W2 = W2.tolil()
    j = 0
    for k in range(n3):
        for i in range(n1):
            idx = k * n12 + j * n1 + i
            W2[idx,:] = W2[idx+n1,:]
    
    j = n2-1
    for k in range(n3):
        for i in range(n1):
            idx = k * n12 + j * n1 + i
            W2[idx,:] = W2[idx-n1,:]
    
    #n3 direction
    W3 = sparse.diags(-2.0*coeff, shape=(n123, n123))
    W3 += sparse.diags(coeff, offsets=n12, shape=(n123, n123))
    W3 += sparse.diags(coeff, offsets=-n12, shape=(n123, n123))
    W3 = W3.tolil()
    k = 0
    for j in range(n2):
        for i in range(n1):
            idx = k * n12 + j * n1 + i
            W3[idx,:] = W3[idx+n12,:]
    
    k = n3-1
    for j in range(n2):
        for i in range(n1):
            idx = k * n12 + j * n1 + i
            W3[idx,:] = W3[idx-n12,:]
    
    #sum up
    W = a1 / d1 / d1 * W1 + a2 /d2 / d2 * W2 + a3 /d3 / d3 * W3
    
    return W.tocsc()

def lcv_lambda(a0,step, num):
    a=np.zeros(num,np.float32)
    for i in range(num):
        a[i]=a0/(step**i)
    return a

def fun(a,b):
    return a+b
