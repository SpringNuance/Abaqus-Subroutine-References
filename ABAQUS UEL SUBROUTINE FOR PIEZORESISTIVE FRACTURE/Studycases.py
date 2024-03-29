############################################################################################
# Leonel Quinteros, Enrique García-Macías, and Emilio Martínez-Pañeda,                     #
# Electromechanical phase-field fracture modelling of piezoresistive CNT-based composites, #
# Computer Methods in Applied Mechanics and Engineering,                                   #
# Volume 407,                                                                              #
# 2023,                                                                                    #
# 115941,                                                                                  #
# ISSN 0045-7825,                                                                          #
# https://doi.org/10.1016/j.cma.2023.115941.                                               #
############################################################################################
#%% Import libraries

import numpy as np


#%% Micromechanical parameters

# Matrix
# Young´s modulus 
EMatrix = 2.5e9  #[MPa]          
# Poisson´s ratio
NuMatrix = 0.28  #[-]              

# MWCNTs
# Isotropic material
# Young´s modulus 
ECnt = 700e9  #[MPa]               
# Poisson´s ratio
NuCnt = 0.3  #[-]                 

# Length 
Lcnt = 3.20995854347668e-6  #[m]               
# Diameter 
Dcnt = 10.3538294770233e-9  #[m]             
# Young´s modulus
EInt = 2.17e9   #[MPa]              
# Poisson´s ratio
NuInt = NuMatrix                  
# Interphase thickness nm
Intert = 31e-9  #[m]                     
# 1- Soft interphase, 2- Hard interphase
type_inter = 1               
# Volume fraction
fc = 0.01 
# Pristine material Critical energy release rate
G0 = 133# [J/m2]
# CNT ultimate strength
SigmaUlt = 35e9  # [Pa]
# Theta angle integration Fracture Energy
ThetaMin = 0  # [rad]
ThetaMax = np.pi / 2  # [rad]
# Snubbing coeff
mu = 0  
# Orientation limit angle
Ac = 0.083 
# CNT interfatial strength [Pa] (epoxy-CNT)
TauInt = 47e6  

# %% Estimate using a MeanField Homogenisation approach
#    Young's modulus and Poisson ratio of the composite

from MfhFunctions import MFH

HMethod = MFH(Lcnt, Dcnt, EMatrix, NuMatrix, ECnt,
                  NuCnt, EInt, NuInt, Intert, fc)

E, nu = HMethod.ComputeMechanicalProps()

# %% Estimate fracture energy

from Fracture_energy import FractureEnergy

GM = FractureEnergy(G0, Lcnt, Dcnt, SigmaUlt, TauInt, Ac, mu,
                        ECnt, ThetaMin, ThetaMax, fc, p=0.5, q=0.5,
                        NInter=100, Ntheta=200)

Gc = GM.EnergyReleaseRate()

# %%  Electrical conductivity and piezoresistive coefficients
# Please refer to github.com/EnriqueGarMac

cond = 0.0337437  # [S/m]

pi11 = 1.2017119765168127

pi12 = 1.2015765068624864

pi44 = 6.773482716315726e-05

# %% Case studies
# Import Abaqus modules
import os
from abaqus import *
from abaqusConstants import *
from caeModules import *
from abaqus import getInput
from odbAccess import *
# Import Study cases
from Studies import *
# PLEASE SELECT THE DIRECTORY
os.chdir(
    r"/rds/general/user/ldq20/home/case_study_length/L_scale_study"
)

Vimp = 10.0  # [v] Applied Voltage

Dimp = 0.0001  # [m] Applied Displacement

Inc = 0.005  # Increment of steps

NotchedAngle = 0 # [-]

k_value = 50  #

n_value = 6  #

meshsize1 = 0.001

meshsize2 = meshsize1*3 

case_number = 1 # Could be from 1 to 5

file_name = 'CaseStudy'

case_study(
        Vimp,
        Dimp,
        Inc,
        E,
        nu,
        Gc,
        NotchedAngle,
        cond,
        pi11,
        pi12,
        pi44,
        k_value,
        n_value,
        file_name,
        case_number,
        meshsize1,
        meshsize2
        )


