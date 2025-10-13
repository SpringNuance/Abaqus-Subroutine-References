#####################################################################################
# THIS CODE CALCULATES THE AVERAGE SURFACE FACTOR
#####################################################################################
"""
@authors: Sara Jimenez Alfaro & Emilio Martinez Paneda
"""
# Last update: 09/05/2025
from math import *
import numpy as np 


def Ks_calculator(a, b, ivector,name_analysis,sigma,error):

    CC, MC, VC, ME, s1e6_dist = [], [], [], [], []

    for sim in ivector:
        name_document = str(name_analysis)+f'_Sim{sim}_numerics.dat'
        iresults = np.loadtxt(name_document)
        iresults = np.vstack(([0,0,0],iresults))
        CC.append(iresults[iresults.shape[0]-1,0])

    for (j,Nc) in enumerate(CC):
        MC.append(np.mean(CC[0:j+1])) # expected value
        VC.append(np.std(CC[0:j+1]))  # standard deviation
    for ci in CC:
        A = sigma*ci**(0.0442)
        s1e6_dist.append(A*(1e6)**(-b)) # stress distribution for 1 million cycles 

    for (j,Nc) in enumerate(CC[29:]):
        k = j + 29
        ME.append(1.96*np.std(s1e6_dist[0:k+1])/sqrt(k+1)) # margin of error
        
    me_global = (1.96*np.std(s1e6_dist)/sqrt(len(s1e6_dist))) # Margin of error of the entire sample (in stress for 1 million)
    s1e6_pos = np.mean(s1e6_dist) + me_global # Stress at 1 million cycles plus the margin of error
    s1e6_neg = np.mean(s1e6_dist) - me_global # Stress at 1 million cycles minus the margin of error
    s1e6     = np.mean(s1e6_dist)             # Stress at 1 million cycles
    me_por = me_global / s1e6 * 100  # percentage of the margin of error with respect to the sample size 
    max_muestra = np.where(me_por < error)[0] # number of cases for which the margin of error is lower than 5%
    
    if len(max_muestra) == 0:
        print('Increase the sample size')
    else:
        max_muestra = max_muestra[0]+30

    surface_factor =  s1e6/(a*1.e6**(-b))
    surface_factor_pos =  s1e6_pos/(a*1.e6**(-b))
    surface_factor_neg =  s1e6_neg/(a*1.e6**(-b))

    return surface_factor, surface_factor_pos, surface_factor_neg, me_por, MC[len(MC)-1]