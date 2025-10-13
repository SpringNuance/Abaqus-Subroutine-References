#####################################################################################
# THIS CODE IS THE INPUT FILE TO ESTIMATE THE FATIGUE LIFE OF THE ROUGH SPECIMEN
#####################################################################################
"""
@authors: Sara Jimenez Alfaro & Emilio Martinez Paneda
"""
# Last update: 08/05/2025
import sys
import numpy as np
sys.path.append("./phasefield_pycodes")
import time
from math import *
from mpi4py import MPI  

# Dolfinx version v0.7.3

#####################################################################################
# Initial parameters 
#####################################################################################
Parameters_data = { "Load":340.,"Gc":18.24,"l0":2.9,"E": 200e3,"nu": 0.3, "model":"AT1","NCycles":5000,'Njump':1,"a":485.9,"b":0.0442, "Ra":1.5, "lcorl":30}       
simulation_vector = np.linspace(start = 1, stop = 1, num = 1)
rmsr = 1.25*Parameters_data.get('Ra')/1000
Parameters_mesh = {"m0":min(Parameters_data.get("lcorl")/1000, Parameters_data.get("l0"))/5, "m1":1.0, "m2":0.05, "n":20, "d":0.8}

#####################################################################################
# Iteration procedure
#####################################################################################
for sim in simulation_vector:

    name_output = f'roughtest_R{Parameters_data.get("Ra")}_L{Parameters_data.get("lcorl")}_Sim{sim}'
    dir_output  = f'roughtest_R{Parameters_data.get("Ra")}_L{Parameters_data.get("lcorl")}/roughtest_R{Parameters_data.get("Ra")}_L{Parameters_data.get("lcorl")}_Sim{sim}'

    from mesh_generation.roughmesh import roughmesh
    mesh, cells, facets, coords_point = roughmesh(rmsr,Parameters_data.get("lcorl")/1000,Parameters_mesh,dir_output)

    #####################################################################################
    # Phase Field iteration (Force control)
    #####################################################################################
    from pfsolver_fatigue import phasefield_solver
    phasefield_solver(name_output,mesh, cells, facets, Parameters_data, coords_point)


#############################################################
## FINAL CALCULATION
#############################################################
from Ks_calculator import Ks_calculator

error = 5
Ks_R, Ksp_R, Ksn_R, me_R = np.zeros((len(sim_R),)), np.zeros((len(simulation_vector),)), np.zeros((len(simulation_vector),)), np.zeros((len(simulation_vector),))
for (j,S) in enumerate(sim_R):
    name_analysis = f'output/roughtest_R{Parameters_data.get("Ra")}_L{Parameters_data.get("lcorl")}'
    Ks_R[j], Ksp_R[j], Ksn_R[j], me_R[j], NF = Ks_calculator(Parameters_data.get("a"), Parameters_data.get("b"), simulation_vector,name_analysis,Parameters_data.get("Load"),error)


# print results for the experimental comparison
print('**************************************')
print('Surface factor')
print('**************************************')
print('Minimum and maximum Ks: ', np.min(Ks_R), np.max(Ks_R))
print('Margin of error for the minimum and maximum Ks: ', me_R[np.argmin(Ks_R)], me_R[np.argmax(Ks_R)])
print('Number of cycles to failure obtained in the simulation:', NF)

    