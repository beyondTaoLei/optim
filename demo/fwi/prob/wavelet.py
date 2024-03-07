#!/usr/bin/env python3

import os
import sys
import obspy
import numpy as np
import pandas as pd
from obspy.io.segy.segy import SEGYTraceHeader, SEGYBinaryFileHeader

def wavelet(shape, nt, delta, tshift, fc, a):
    amp = np.zeros(nt,  np.float32)
    t=np.arange(1, nt+1)*delta
    PI = np.pi
    ts=1.0/fc
    
    if shape == 1: # New Ricker Wavelet, equal to SOFI2D
        tau = PI*(t-1.5*ts-tshift)/ts
        tau2 = tau*tau
        amp[:] = (1.0-2.0*tau2)*np.exp(-tau2)
    elif shape == 2:
        amp[:] = np.sin(2.0*PI*(t-tshift)*fc) - 0.5*np.sin(4.0*PI*(t-tshift)*fc)
        amp[t<tshift] = 0.0
        amp[t>(tshift+ts)] = 0.0
    elif shape == 3: # sinus raised to the power of three
        amp[:] = np.power(np.sin(PI*(t+tshift)/ts),3.0)
        amp[t<tshift] = 0.0
        amp[t>(tshift+ts)] = 0.0
    elif shape == 4: # first derivative of a Gaussian
        ts=1.2/fc
        ag  = PI*PI*fc*fc
        tts = t-ts
        amp[:] = -2.0 * ag * tts * np.exp(-ag*tts*tts)
    elif shape == 5: # Bandlimited Spike
        it = round(tshift/delta)
        amp[it] = 1.0    
    elif shape == 6: #integral of sinus raised to the power of three
        cos1 = np.cos(PI*(t-tshift)/ts)
        amp[:] = (ts/(0.75*PI))*(0.5-0.75*cos1+0.25*np.power(cos1,3.0))
        amp[t<tshift] = 0.0
        amp[t>(tshift+ts)] = ts/(0.75*PI)
    elif shape == 7: #source wavelet from file SOURCE_FILE
        pass
    else:
        raise ValueError('Which source-wavelet ? ')
    
    return amp*a

if __name__ == "__main__":
    #Input paras
    fevn_info   =sys.argv[1] #list/eventinfo.csv
    shape       =int(sys.argv[2]) #1
    nt          =int(sys.argv[3]) #6500
    delta       =float(sys.argv[4]) #0.001
    fout        =sys.argv[5] #source
    
    # interpretation
    df = pd.read_csv(fevn_info)
    name    = list(df['station'])
    X       = list(df['x'])
    Y       = list(df['y'])
    tshift  = list(df['tshift'])
    fc      = list(df['fc'])
    amp0    = list(df['amp'])
    nsrc    = len(name)
    
    # generate the wavelets, and save to su-format file
    for i in range(nsrc):
        print('wavelet: %4d of %4d'%(i+1, nsrc))
        scalel = 1000.0
        scalco = 1000.0
        itr = 1
        wav = wavelet(shape, nt, delta, tshift[i], fc[i], amp0[i])
        # create Stream
        stream = obspy.Stream()
        trace = obspy.Trace(data=wav)
        # Attributes in trace.stats will overwrite everything in
        # trace.stats.segy.trace_header
        trace.stats.starttime = obspy.UTCDateTime(1970,1,1,00,00,00)
        trace.stats.delta = delta
        
        # format: obspy.io.segy.header.TRACE_HEADER_FORMAT
        trace.stats.su = {}
        trace.stats.su.trace_header = SEGYTraceHeader()
        tr_header = trace.stats.su.trace_header
        tr_header.trace_sequence_number_within_line = itr
        tr_header.trace_identification_code = 1 # trace identification code: 1=seismic
        tr_header.receiver_group_elevation = int(Y[i] * scalel) # not use
        tr_header.water_depth_at_source = 1000 # not use
        tr_header.scalar_to_be_applied_to_all_elevations_and_depths = - int(scalel)
        tr_header.scalar_to_be_applied_to_all_coordinates = -int(scalco)
        tr_header.source_coordinate_x   = int(X[i] * scalco) #Sx
        tr_header.source_coordinate_y   = int(Y[i] * scalco) #Sy
        tr_header.group_coordinate_x    = int(X[i] * scalco) #Rx=Sx
        tr_header.group_coordinate_y    = int(Y[i] * scalco) #Ry=Sy
        tr_header.number_of_samples_in_this_trace = nt 
        tr_header.sample_interval_in_ms_for_this_trace = int(delta * 1.0e6) # sample interval in micro-seconds
        
        # Add trace to stream
        stream.append(trace)
        stream.stats = obspy.core.AttribDict()
        stream.stats.textual_file_header = 'Textual Header!'
        stream.stats.binary_file_header = SEGYBinaryFileHeader()
        stream.stats.binary_file_header.trace_sorting_code = 5
        fnm = os.path.join(fout, name[i]+'.su')
        stream.write(fnm, format='SU', byteorder='<')
        
print("Finished...", __file__)