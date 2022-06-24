#!/usr/bin/env python3
"""
Usage:
    python *.py
"""
import os
import sys
import obspy
import pandas as pd
import obspy.signal.filter as flt

#Input paras
fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
fsource     =optim['fsource']
fevn_info   =optim['fevn_info']
fc          =optim['fc']   #low cut-off frequency
order       =optim['order']#4
zero_phase  =optim['zero_phase']
zero_phase =False if zero_phase==0 else True

events=pd.read_csv(fevn_info)['station']
for i in range(len(events)):
    fin =os.path.join(fsource, 'wavelet_shot'+str(i+1)+'.su')
    fout =os.path.join(optim['fdio'], events[i],'wavelet.su')
    if os.path.exists(fin):
        print(fin+" Butterworth filtering( FC ORDER= %.2f %d )"%(FC, ORDER))
        wav=obspy.read(fin)[0]
        delta=wav.stats.delta
        #Apply a zero-phase lowpass to the trace
        wav.data[:]=flt.lowpass(wav.data,FC,1.0/delta, ORDER, zerophase=zero_phase)
        wav.write(fout)
    else:
        raise ValueError("no file: "+fin)

print("Finished...", __file__)
        
