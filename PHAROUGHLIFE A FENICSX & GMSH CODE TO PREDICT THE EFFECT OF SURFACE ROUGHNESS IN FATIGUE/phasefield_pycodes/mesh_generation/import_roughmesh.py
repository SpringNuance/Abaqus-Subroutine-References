#####################################################################################
# THIS CODE IMPORT THE CREATED ROUGH PROFILE AND DEFINE THE CORRESPONDING PHYSICAL GROUPS
#####################################################################################
"""
@authors: Sara Jimenez Alfaro & Emilio Martinez Paneda
"""
# Last update: 08/05/2025
import gmsh
import dolfinx
from math import*
from mpi4py import MPI
import numpy as np
import math

def import_mesh(mesh_output, Parameters):

    model_rank = 0
    gmsh.initialize()
    gdim = 2
    mesh_comm = MPI.COMM_WORLD

    coords_point = []
    if mesh_comm.rank == 0:

        #####################################################################################
        # GEOMETRICAL AND MESH PARAMETERS
        #####################################################################################
        geomlines = 15 # Number of geometrical boundaries 
        m0,m1,m2 = Parameters.get("m0"), Parameters.get("m1"), Parameters.get("m2")
        a_ellipse = sqrt(12.2**2/(1-((6.9-3.8-25.9)/25.9)**2))
        n, d = Parameters.get("n"), Parameters.get("d")

        # Directory of the rough mesh
        ruta_malla2 = mesh_output+"/mesh_rough.ply2"
        # Extract the vertex and elements in the example mesh
        vertex, face = read_ply2(ruta_malla2)
        
        tagi = 1
        # Additional inner lines
        tags_l18, ptags_l18 = [], [] 
        tags_l16, ptags_l16 = [], [] 
        tags_l17, ptags_l17 = [], [] 
        # Physical groups for the boundary conditions
        ptags_left,tags_left,tags_right   = [],[], []
        
        x,y,z = vertex[0,0],vertex[0,1],vertex[0,2]
        mi = m1
        gmsh.model.geo.addPoint(x,y,z,mi,tagi)
        tags_left.append(tagi)
        for i in range(1,geomlines+1): 
            nodes = np.loadtxt(mesh_output+"/boundnodes/nodes_boundary_l{}.dat".format(i), dtype=int)
            for (j,nodo_bound) in enumerate(nodes[1:len(nodes)]):
                # If we are in the last line we need to join with the first point
                if i == geomlines and j==len(nodes)-2:
                    tagi = tagi + 1
                    gmsh.model.geo.addLine(1,tagi-1)
                else:
                    tagi = tagi + 1
                    # Refinement
                    if i == 3 or i == 11:
                        mi = m0
                    # If any of the lines is inside the physical group save it
                    if i == 7:
                        tags_right.append(tagi-1)
                    if i == 14 or i==15:
                        tags_left.append(tagi-1)
                    if i == 15:
                        if j == 0:
                            ptags_left.append(tagi)
                    # If any point belongs to the inner lines, we store it.
                    if i == 2 or i == 11:
                        if j == len(nodes)-2:
                            mi = m2 # Refinement in these places
                            ptags_l17.append(tagi)
                            tags_l17.append(tagi-1)
                    if i == 3 or i == 10:
                        if j == len(nodes)-2:
                            mi = m0 # Refinement in these places
                            ptags_l16.append(tagi)
                            tags_l16.append(tagi-1)
                    if i == 4 or i == 9:
                        if j == len(nodes)-2:
                            mi = m0 # Refinement in these places
                            ptags_l18.append(tagi)
                            tags_l18.append(tagi-1)
                    # Coordinates of the point to save
                    x,y,z = vertex[nodo_bound,0],vertex[nodo_bound,1],vertex[nodo_bound,2]
                    gmsh.model.geo.addPoint(x,y,z,mi,tagi)
                    # It is the line segment that connects two consecutive points belonging to the inner lines
                    gmsh.model.geo.addLine(tagi,tagi-1,tagi-1)
                mi = m1

        p16 = gmsh.model.geo.addPoint(32.5-a_ellipse/n,6.9+3.8+25.9-25.9*sqrt(1-(a_ellipse/n)**2/a_ellipse**2)-d,0,m2,tagi)
        p17 = gmsh.model.geo.addPoint(32.5-a_ellipse/n,6.9-3.8-25.9+25.9*sqrt(1-(a_ellipse/n)**2/a_ellipse**2)+d,0,m2,tagi+1)
        p18 = gmsh.model.geo.addPoint(32.5+a_ellipse/n,6.9+3.8+25.9-25.9*sqrt(1-(a_ellipse/n)**2/a_ellipse**2)-d,0,m2,tagi+2)
        p19 = gmsh.model.geo.addPoint(32.5+a_ellipse/n,6.9-3.8-25.9+25.9*sqrt(1-(a_ellipse/n)**2/a_ellipse**2)+d,0,m2,tagi+3)
        p21 = gmsh.model.geo.addPoint(32.5,6.9-3.8+d,0,m2)
        p20 = gmsh.model.geo.addPoint(32.5,6.9+3.8-d,0,m2)
        p23 = gmsh.model.geo.addPoint(32.5,6.9,0,m2)
        # Include the interior lines
        gmsh.model.geo.addLine(ptags_l16[0],p21,tagi)
        gmsh.model.geo.addLine(ptags_l17[0],p17,tagi+1)
        gmsh.model.geo.addLine(ptags_l18[0],p19,tagi+2)
        l19 = gmsh.model.geo.addLine(p17,p16,tagi+3)
        l20 = gmsh.model.geo.addLine(p16,ptags_l17[1],tagi+4)
        l21 = gmsh.model.geo.addLine(p21,p23,tagi+5)
        l22 = gmsh.model.geo.addLine(p20,ptags_l16[1],tagi+6)
        l23 = gmsh.model.geo.addLine(p19,p18,tagi+7)
        l24 = gmsh.model.geo.addLine(p18,ptags_l18[1],tagi+8)
        l25 = gmsh.model.geo.addLine(p16,p20,tagi+9)
        l26 = gmsh.model.geo.addLine(p17,p21,tagi+10)
        l27 = gmsh.model.geo.addLine(p21,p19,tagi+11)
        l28 = gmsh.model.geo.addLine(p20,p18,tagi+12)
        l29 = gmsh.model.geo.addLine(p23,p20,tagi+13)
        # Generate a curve loop
        cloop1 = gmsh.model.geo.addCurveLoop(np.concatenate((np.arange(-1,-tags_l17[0]-1, -1), [tagi+1,tagi+3,tagi+4], np.arange(-tags_l17[1]-1,-tagi, -1))).tolist())
        surface1 = gmsh.model.geo.addPlaneSurface([cloop1],1)
        cloop2 = gmsh.model.geo.addCurveLoop(np.concatenate((np.arange(-tags_l17[0]-1,-tags_l16[0]-1, -1), [tagi,-(tagi+10),-(tagi+1)])).tolist())
        surface2 = gmsh.model.geo.addPlaneSurface([cloop2],2)
        cloop3 = gmsh.model.geo.addCurveLoop(np.concatenate((np.arange(-tags_l16[0]-1,-tags_l18[0]-1, -1), [tagi+2,-(tagi+11),-tagi])).tolist())
        surface3 = gmsh.model.geo.addPlaneSurface([cloop3],3)
        cloop4 = gmsh.model.geo.addCurveLoop(np.concatenate((np.arange(-tags_l18[0]-1,-tags_l18[1]-1, -1), [-(tagi+8),-(tagi+7),-(tagi+2)])).tolist())
        surface4 = gmsh.model.geo.addPlaneSurface([cloop4],4)
        cloop5 = gmsh.model.geo.addCurveLoop([l26,l21,l29,-l25,-l19])
        surface5 = gmsh.model.geo.addPlaneSurface([cloop5],5)
        cloop6 = gmsh.model.geo.addCurveLoop([l27,l23,-l28,-l29,-l21])
        surface6 = gmsh.model.geo.addPlaneSurface([cloop6],6)
        cloop7 = gmsh.model.geo.addCurveLoop(np.concatenate((np.arange(-tags_l18[1]-1,-tags_l16[1]-1, -1), [-(tagi+6),tagi+12,tagi+8])).tolist())
        surface7 = gmsh.model.geo.addPlaneSurface([cloop7],7)
        cloop8 = gmsh.model.geo.addCurveLoop(np.concatenate((np.arange(-tags_l16[1]-1,-tags_l17[1]-1, -1), [-(tagi+4),tagi+9,tagi+6])).tolist())
        surface8 = gmsh.model.geo.addPlaneSurface([cloop8],8)
        

        # Generate the mesh for the original half of the geometry.
        gmsh.model.geo.synchronize()

        # Definition of physical groups 3D mesh
        # Domain
        gmsh.model.addPhysicalGroup(gdim, [surface1,surface2,surface3,surface4,surface5,surface6,surface7,surface8],1) 
        gmsh.model.setPhysicalName(gdim, 1, "domain")
        pl2=gmsh.model.addPhysicalGroup(1, tags_left,3)
        gmsh.model.setPhysicalName(1, pl2, "left")
        pl2=gmsh.model.addPhysicalGroup(1, tags_right,4)
        gmsh.model.setPhysicalName(1, pl2, "right")
            
        # Generation of a 2D mesh 
        gmsh.model.mesh.generate(gdim)

        # Visualizar la malla 
        # gmsh.fltk.run()

    
    #####################################################################################
    # MESH EXPORTATION
    #####################################################################################
    mesh, mt, ft = dolfinx.io.gmshio.model_to_mesh(gmsh.model,mesh_comm,model_rank,gdim)
        

    def get_node_coordinates(node_tags):
        coords = []
        for node_tag in node_tags:
            x, y, z = gmsh.model.mesh.getNode(node_tag)[0:3]
            coords.append((x, y, z))
        return coords
    
    coords_point = []

    if mesh_comm.rank == 0:
        coords_point = get_node_coordinates(ptags_left)[0][0]
        for dest in range(1, mesh_comm.size):
            mesh_comm.send(coords_point, dest=dest)
    else:
        # Receive the value in the other processes.
        coords_point = mesh_comm.recv(source=0)

    
    # Close GMSH
    gmsh.finalize()

    return mesh,mt,ft, coords_point


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

    
    