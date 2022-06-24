#!/usr/bin/env python3
"""
Usage:
    python *.py
"""
import os
import sys
import json
import obspy
from misfit_func import misf_waveform

#Input paras
flag    =int(sys.argv[1])#1: output adj_src; 0: no adj_src
lnorm   =int(sys.argv[2])#1: L2 norm
fsyn    ="syn_final.su"
fobs    ="obs_final.su"
fmisfit ="misfit.json"
fadj    ="adj_shot1.su"

if not (os.path.exists(fsyn) and os.path.exists(fobs)):
    raise ValueError("no file: "+fsyn +"\n or "+fobs)
syn=obspy.read(fsyn)
obs=obspy.read(fobs)
if lnorm==1:
    misfit, energy, adj=misf_waveform(syn, obs)
else:
    raise ValueError("we will implement other norm in future!")

if flag==1:
    adj.write(fadj)
misf_paras={
    'misfit'    :misfit,
    'energy'    :energy
}
fout = open(fmisfit,'w')
fout.write(json.dumps(misf_paras,indent=4))
fout.close()

print("Finished...", __file__)
        
