#!/usr/bin/env python3
import numpy as np
            
def angle_restart(grad, desc, tau=-0.02):
    """
    @article{modrak_seismic_2016,
            title = {Seismic waveform inversion best practices: regional, global and exploration test cases},
            volume = {206},
            issn = {0956-540X, 1365-246X},
            url = {https://academic.oup.com/gji/article-lookup/doi/10.1093/gji/ggw202},
            doi = {10.1093/gji/ggw202},
            shorttitle = {Seismic waveform inversion best practices},
            pages = {1864--1889},
            number = {3},
            journaltitle = {Geophysical Journal International},
            author = {Modrak, Ryan and Tromp, Jeroen},
            urldate = {2019-01-04},
            date = {2016-09-01},
            langid = {english},
            keywords = {global seismology}
    }
    """
    #ANGLE RESTART CONDITION
    flag=0
    q0 =np.inner(grad, desc)
    ng =np.linalg.norm(grad)
    nd =np.linalg.norm(desc)
    if q0/ng/nd>tau:
        flag=1
        print('RESTART LINESEARCH BASED ON ANGLE RESTART CONDITION')
        print('q0 :',q0)
        print('ng :',ng)
        print('nd :',nd)
        print('q0/ng/nd,tau :',q0/ng/nd,tau)
    return flag
    
        