#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 15:28:44 2019
------------------------------------------------
                 PRUEBA DE SNR
------------------------------------------------
@author: soporte
"""
import numpy as np
import drawFigures as df
import os
 
x = np.linspace(0.001,50,5000)
SNR = np.zeros(len(x))
STD = np.zeros(len(x))
k = 0

for i in x:
    f = np.array([1]*500) # Function
    noise  = np.random.normal(0,i,f.shape) + 1j*np.random.normal(0,i,f.shape) # Noise function 
    s = f+noise

    # Getting SNR value
    SNR[k] = 1/(noise.std()**2)

    # Getting STD of the phase 
    phase = np.angle(s)
    STD[k] = phase.std()
    
    k += 1

SNR = 10*np.log10(SNR) # Logarithmic scale 

drc = os.getcwd()+"/Results/Otros/snr_test.png"
df.simple_Uplot(SNR,STD,title=r"$\sigma_\phi \ vs \ SNR$",axis_X=r"$SNR(dB)$",axis_Y=r"$\sigma_\phi$",direction=drc)

