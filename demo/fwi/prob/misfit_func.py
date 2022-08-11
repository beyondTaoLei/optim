#!/usr/bin/env python3
import numpy as np

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

def misf_waveform(syn, obs):
    """
    L2 norm misfit
    misfit = 0.5<Ru-d_obs|Ru-d_obs>
    adj = Ru-d_obs
    """
    ntr =syn.count()
    delta =syn[0].stats.delta
    adj=obs.copy()
    misfit, energy =0., 0.
    for itr in range(ntr):
        adj[itr].data[:]=syn[itr].data-obs[itr].data
        misfit +=np.inner(adj[itr].data, adj[itr].data)
        energy +=np.inner(obs[itr].data, obs[itr].data)
    misfit = 0.5*misfit*delta
    energy = 0.5*energy*delta
    
    return misfit, energy, adj
