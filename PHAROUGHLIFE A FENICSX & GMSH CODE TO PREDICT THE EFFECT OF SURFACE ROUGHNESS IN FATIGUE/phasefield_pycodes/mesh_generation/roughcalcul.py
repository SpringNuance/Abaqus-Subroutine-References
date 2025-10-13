
#####################################################################################
# Generate unstructured mesh of general two-dimensional rough surface
#####################################################################################
"""
@authors: Sara Jimenez Alfaro & Emilio Martinez Paneda
"""
# Last update: 08/05/2025
import numpy as np
from scipy.spatial.distance import pdist, squareform
import os
import logging
import sys
import math
from mpi4py import MPI

def roughcalcul(fname, corl, rmsr, mesh_output, mesh_comm):
    
    if mesh_comm.rank == 0:
        
        # File stem and extension
        stem, ext = os.path.splitext(os.path.basename(fname))

        # Read the nodes of the boundary that were saved in the smoothmesh.py
        dat_filename = mesh_output+"/boundnodes/nodes_boundary.dat"
        # Load the labels of the nodes in a numpy vector
        nodes = np.sort(np.loadtxt(dat_filename, dtype=int))
        
        # Extract the vertex and elements in the example mesh
        vertex, face = read_ply2(fname)   
        vertex_nodes = vertex[nodes]
        
        # Compute Eulidean distance matrix
        distmat = squareform(pdist(vertex_nodes))
        # Compute the matrix R
        Rmat = np.exp(-0.5*(distmat/corl)**2)
        # Cholesky decomposition 
        Lmat = chol_eigh(Rmat)
        # Generate a white noise
        gaussian_white_noise = np.random.randn(vertex_nodes.shape[0])
        # Compute the perpendicular distance in the roughness
        height = rmsr * (Lmat @ gaussian_white_noise)
        # Get the normal vectors in the mesh
        normal = get_normal(vertex, face, nodes)
        # Compute the changes in the vertex
        vertex_rough = vertex_nodes + (height * normal.T).T
        # Change the global coordinates
        for (i,j) in enumerate(nodes):
            vertex[j,:] = vertex_rough[i,:]

        # output filename
        fname_rough = stem + '_rough' + ext
        full_fname_rough = os.path.join(os.path.dirname(fname),fname_rough)
        
        write_ply2(full_fname_rough, vertex, face)

def read_ply2(fname):
    # This function reads the .ply2
    with open(fname, 'r') as f:

        vertex_count = np.fromstring(f.readline(), dtype=int, sep=' ')[0]
        vertex = np.zeros((vertex_count,3), dtype=float)

        face_count = np.fromstring(f.readline(), dtype=int, sep=' ')[0]
        face = np.zeros((face_count,3), dtype=int)

        for i in range(vertex_count):
            vertex[i,:] = np.fromstring(f.readline(), dtype=float, sep=' ')

        for i in range(face_count):
            face[i,:] = np.fromstring(f.readline(),dtype=int,sep=' ')[1:]

    return vertex, face

def write_ply2(fname, vertex, face):
    # This function writes the .ply2
    with open(fname, 'w') as f:
        vertex_count = vertex.shape[0]
        face_count = face.shape[0]
        f.write('{}\n{}\n'.format(vertex_count, face_count))
        for i in range(vertex_count):
            f.write('{} {} {}\n'.format(vertex[i,0], vertex[i,1], vertex[i,2]))
        for i in range(face_count):
            f.write('3 {} {} {}\n'.format(face[i,0], face[i,1],face[i,2]))

def get_normal(vertex, face, nodes):

    """Return normal vectors in a 2d mesh"""

    # Extract the elements where there are more than two nodes in the boundary (boundary elements)
    elements_boundary = []
    
    # Read the matrix row by row
    for row in face:
        # Determine the number of elements in the row that are in the vector "nodes"
        boundnodes_in_element = sum(1 for elem in row if elem in nodes)
        # If there is more than one element of the row in "nodes", save that row
        if boundnodes_in_element > 1:
            elements_boundary.append(row)
    
    elements_boundary = np.vstack(elements_boundary)
    
    # Build the matriz with the normal vectors outward-facing normals to the contour (two per node).
    normal_elements_boundary = np.zeros((len(nodes),6))

    for i in range(len(elements_boundary)):

        # Extract the two nodes inside the contour and the node outside
        nodbound = elements_boundary[i,:][np.isin(elements_boundary[i,:], nodes)]
        nopbound = elements_boundary[i,:][~np.isin(elements_boundary[i,:], nodes)]      
        
        # Coordinates of the nodes that are in the boundary
        coord1, coord2 = vertex[nodbound[0]],vertex[nodbound[1]]
        # Coordinates of the nodes that are NOT in the boundary
        coord3 =  vertex[nopbound[0]]   

        # OUTWARD NORMAL VECTOR CENTERED IN THE POINT 1
        # Case: Vertical line
        if coord2[0]-coord1[0] == 0:
            # No rotation is needed. Y-coordinate of point 3 in the rotated coordinate system.
            x3p = (coord3[0]-coord1[0])
            # Normal vector in the translated coordinate axis.
            vp  = np.array([-x3p,0.0,0.0])
            # Perpendicular, outward, and unit vector centered at point 1.
            vp  = vp / np.linalg.norm(vp)
        else:
            # Slope between points 1 and 2
            mp  = (coord2[1]-coord1[1])/(coord2[0]-coord1[0])
            # Rotation matrix of the system of coordinates
            Lth  = np.array([[math.cos(math.atan(mp)),-math.sin(math.atan(mp)),0.0],[math.sin(math.atan(mp)),math.cos(math.atan(mp)),0.0],[0.0,0.0,1.0]])
            # Y-coordinate of point 3 in the rotated coordinate system
            y3pp = math.sin(math.atan(mp)) * (coord3[0]-coord1[0]) + math.cos(math.atan(mp)) * (coord3[1]-coord1[1])
            # Normal vector in the rotated coordinate system
            vp  = np.dot(np.transpose(Lth),np.array([0.0,-y3pp,0.0]))
            # Perpendicular vector, outward and unitary, centered in the point 1
            vp  = vp / np.linalg.norm(vp) 
        
        # Save the normal in a global matrix
        if np.all(normal_elements_boundary[np.where(nodes == nodbound[0])[0],:] == 0):
            normal_elements_boundary[np.where(nodes == nodbound[0])[0],0:3] = vp
        else:
            normal_elements_boundary[np.where(nodes == nodbound[0])[0],3:] = vp

        # OUTWARD NORMAL VECTOR CENTERED IN THE POINT 2
        # Case: Vertical line
        if coord1[0]-coord2[0] == 0:
            # No rotation is needed. Y-coordinate of point 3 in the rotated coordinate system.
            x3p = (coord3[0]-coord2[0])
            # Normal vector in the translated coordinate axis.
            vp  = np.array([-x3p,0.0,0.0])
            # Perpendicular, outward, and unit vector centered at point 1.
            vp  = vp / np.linalg.norm(vp)
        else:
            # Slope between points 1 and 2
            mp  = (coord1[1]-coord2[1])/(coord1[0]-coord2[0])
            # Rotation matrix of the system of coordinates
            Lth  = np.array([[math.cos(math.atan(mp)),-math.sin(math.atan(mp)),0.0],[math.sin(math.atan(mp)),math.cos(math.atan(mp)),0.0],[0.0,0.0,1.0]])
            # Y-coordinate of point 3 in the rotated coordinate system
            y3pp = math.sin(math.atan(mp)) * (coord3[0]-coord2[0]) + math.cos(math.atan(mp)) * (coord3[1]-coord2[1])
            # Normal vector in the rotated coordinate system
            vp  = np.dot(np.transpose(Lth),np.array([0.0,-y3pp,0.0]))
            # Perpendicular vector, outward and unitary, centered in the point 1
            vp  = vp / np.linalg.norm(vp) 

        # Vector centered in the global coordinates system (for point 2)
        if np.all(normal_elements_boundary[np.where(nodes == nodbound[1])[0],:] == 0):
            normal_elements_boundary[np.where(nodes == nodbound[1])[0],0:3] = vp
        else:
            normal_elements_boundary[np.where(nodes == nodbound[1])[0],3:] = vp

    normat = np.zeros((len(nodes),3))
        
    for i in range(len(normal_elements_boundary)):
        # Obtain the vectors of row i
        v1 = normal_elements_boundary[i, :3]
        v2 = normal_elements_boundary[i, 3:]
        
        # Calculate the intermediate vector (middle)
        v_medio = (v1 + v2) / 2
        
        # Normalize the intermediate vector to obtain a unitary vector
        normat[i,:] = v_medio / np.linalg.norm(v_medio)

    return normat

def chol_eigh(A):

    # Compute Cholesky Decomposition (L) of matrix A. 
    try:
        logging.info('Compute Cholesky decomposition')
        L = np.linalg.cholesky(A)
    except np.linalg.LinAlgError as e:
        logging.info(e)
        logging.info('Compute eigendecomposition')
        eigenvalues, eigenvectors = np.linalg.eigh(A)

        logging.info('Number of negative eigenvalues: {}'.format(
            eigenvalues[eigenvalues < 0.0].shape[0]))
        logging.info('Smallest eigenvalue: {}'.format(np.amin(eigenvalues)))

        eigenvalues[eigenvalues < 0.0] = 0.0
        L = eigenvectors @ np.diag(np.sqrt(eigenvalues))
    
    return L




