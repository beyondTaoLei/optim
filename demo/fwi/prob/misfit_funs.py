#!/usr/bin/env python3
"""
useful functions
Usage:
    [from util import fun] in another *.py file
https://stackoverflow.com/questions/20309456/call-a-function-from-another-file
"""
import numpy as np
from scipy.signal import correlate, correlation_lags
from obspy.signal.differentiate_and_integrate import integrate_cumtrapz

"""
calculates the misfit or adjoint sources according 
to the related algorithm.
INPUT: 
    observed data
    synthetic data
OUTPUT:
    the misfit between two seismograms
    the adjoint source at receiver position
            
"""

# def misf_waveform(syn, obs, delta):
#   """
#   L2 norm misfit
#   misfit = 0.5<Ru-d_obs|Ru-d_obs>
#   adj = Ru-d_obs
#   """
#   # print(syn.shape)
#   ntr, nt = syn.shape
#   adj = np.zeros_like(syn)
#   misfit, energy =0., 0.
#   for itr in range(ntr):
#     adj[itr,:]=syn[itr,:] - obs[itr,:]
#     misfit += np.inner(adj[itr,:], adj[itr,:])
#     energy += np.inner(obs[itr,:], obs[itr,:])
#   # sum up
#   misfit = 0.5*misfit*delta
#   energy = 0.5*energy*delta
#   return misfit, energy, adj

def misf_wfm(syn, obs, delta):
  """
  L2 norm misfit
  misfit = 0.5<Ru-d_obs|Ru-d_obs>
  adj = Ru-d_obs
  """
  adj = syn - obs
  misfit = 0.5 * np.inner(adj, adj) * delta
  energy = 0.5 * np.inner(obs, obs) * delta
  return misfit, energy, adj

def misf_tt(syn, obs, delta, tdiff=None):
  """
  syn and obs : velocity components tapered with specified time window
  traveltime misfit
  misfit = 0.5 (\Delta T)^2
  \Delta T = Tsyn - Tobs
  adj = \frac{\Delta T}{N_i} * \partial_t Dsyn
  N_i = \int si *  \partial^2_t Dsyn dt
  Luo, Y., R. Modrak, and J. Tromp, 2013, Strategies in Adjoint Tomography.
  """
  nt = len(syn)
  if tdiff is None:
    Dsyn = integrate_cumtrapz(syn, delta)
    Dobs = integrate_cumtrapz(obs, delta)
    cc = correlate(Dsyn, Dobs) # T_diff = T_first - T_second
    lags = delta * correlation_lags(nt, nt)
    tdiff = lags[np.argmax(cc)]
  Ni = np.inner(syn, syn) * delta
  adj = tdiff / Ni * syn
  # print(tdiff, Ni)
  misfit = 0.5 * tdiff * tdiff * delta
  return misfit, adj

def fun(a,b):
    return a+b
