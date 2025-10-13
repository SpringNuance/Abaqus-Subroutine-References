#####################################################################################
# THIS CODE COORDINATES THE GENERATION OF THE ROUGH PROFILE
#####################################################################################
"""
@authors: Sara Jimenez Alfaro & Emilio Martinez Paneda
"""
# Last update: 08/05/2025
import sys
import numpy as np
from mesh_generation.roughcalcul import*
from mesh_generation.smoothmesh import smoothmesh
from mesh_generation.import_roughmesh import import_mesh
import os
import shutil
from mpi4py import MPI

def roughmesh(rmsr,corl,Parameters_mesh,nameout):

    mesh_comm = MPI.COMM_WORLD

    #####################################################################################
    # Save the solution
    #####################################################################################
    mesh_output = "./output/"+nameout+"/mesh_output_R{}_l{}".format(rmsr, int(corl))
    if mesh_comm.rank == 0:
        if not os.path.exists(mesh_output):
            os.makedirs(mesh_output,exist_ok=True)
            if not os.path.exists(mesh_output+"/boundnodes"):
                os.makedirs(mesh_output+"/boundnodes",exist_ok=True)
    fname = mesh_output+"/mesh.ply2"
    
    #####################################################################################
    # Smooth mesh
    #####################################################################################
    smoothmesh(Parameters_mesh, mesh_output, mesh_comm)

    #####################################################################################
    # Roughmesh generation
    #####################################################################################
    roughcalcul(fname, corl, rmsr, mesh_output, mesh_comm)
    mesh, cells, facets, coords_point = import_mesh(mesh_output, Parameters_mesh)

    return mesh, cells, facets, coords_point