#####################################################################################
# THIS CODE BUILTS THE NOMINAL SURFACE OF THE ROUGH PROFILE
#####################################################################################
"""
@authors: Sara Jimenez Alfaro & Emilio Martinez Paneda
"""
# Last update: 08/05/2025
# ................ Libraries required
import gmsh
import dolfinx
from math import*
from mpi4py import MPI

# ATTENTION: The goal is to generate a mesh and output the mesh's PLY2 file, a list with the 
# positions where the boundary nodes are located in the PLY2 vertex matrix (in a .dat file)
# and a set of .dat files where the positions in the PLY2 vertex matrix are stored for each 
# set of points corresponding to each boundary line. In this particular case, we have two
# construction points that change the label order.


def smoothmesh(Parameters, mesh_output, mesh_comm):

    model_rank = 0
    gmsh.initialize()
    gdim = 2

    if mesh_comm.rank == model_rank:

        #####################################################################################
        # GEOMETRICAL AND MESH PARAMETERS
        #####################################################################################
        m0,m1,m2 = Parameters.get("m0"), Parameters.get("m1"), Parameters.get("m2")
        a_ellipse = sqrt(12.2**2/(1-((6.9-3.8-25.9)/25.9)**2))
        n, d = Parameters.get("n"), Parameters.get("d")
        
        #####################################################################################
        # GENERATE THE POINTS
        #####################################################################################
        p1 = gmsh.model.geo.addPoint(0,0,0,m1)
        p2 = gmsh.model.geo.addPoint(32.5-12.2,0,0,m1)
        p3 = gmsh.model.geo.addPoint(32.5-a_ellipse/n,6.9-3.8-25.9+25.9*sqrt(1-(a_ellipse/n)**2/a_ellipse**2),0,m0)
        p4 = gmsh.model.geo.addPoint(32.5,6.9-3.8,0,m0)
        p5 = gmsh.model.geo.addPoint(32.5+a_ellipse/n,6.9-3.8-25.9+25.9*sqrt(1-(a_ellipse/n)**2/a_ellipse**2),0,m0)
        p6 = gmsh.model.geo.addPoint(32.5+12.2,0,0,m1)
        p7 = gmsh.model.geo.addPoint(65,0,0,m1)
        p8 = gmsh.model.geo.addPoint(65,2*6.9,0,m1)
        p9 = gmsh.model.geo.addPoint(32.5+12.2,2*6.9,0,m1)
        p10 = gmsh.model.geo.addPoint(32.5+a_ellipse/n,6.9+3.8+25.9-25.9*sqrt(1-(a_ellipse/n)**2/a_ellipse**2),0,m0)
        p11 = gmsh.model.geo.addPoint(32.5,6.9+3.8,0,m0)
        p12 = gmsh.model.geo.addPoint(32.5-a_ellipse/n,6.9+3.8+25.9-25.9*sqrt(1-(a_ellipse/n)**2/a_ellipse**2),0,m0)
        p13 = gmsh.model.geo.addPoint(32.5-12.2,2*6.9,0,m1)
        p14 = gmsh.model.geo.addPoint(0,2*6.9,0,m1)
        p15 = gmsh.model.geo.addPoint(0,6.9,0,m1)
        c1 = gmsh.model.geo.addPoint(32.5,6.9-3.8-25.9,0,m1)
        c2 = gmsh.model.geo.addPoint(32.5,6.9+3.8+25.9,0,m1)
        p16 = gmsh.model.geo.addPoint(32.5-a_ellipse/n,6.9+3.8+25.9-25.9*sqrt(1-(a_ellipse/n)**2/a_ellipse**2)-d,0,m2)
        p17 = gmsh.model.geo.addPoint(32.5-a_ellipse/n,6.9-3.8-25.9+25.9*sqrt(1-(a_ellipse/n)**2/a_ellipse**2)+d,0,m2)
        p18 = gmsh.model.geo.addPoint(32.5+a_ellipse/n,6.9+3.8+25.9-25.9*sqrt(1-(a_ellipse/n)**2/a_ellipse**2)-d,0,m2)
        p19 = gmsh.model.geo.addPoint(32.5+a_ellipse/n,6.9-3.8-25.9+25.9*sqrt(1-(a_ellipse/n)**2/a_ellipse**2)+d,0,m2)
        p21 = gmsh.model.geo.addPoint(32.5,6.9-3.8+d,0,m2)
        p20 = gmsh.model.geo.addPoint(32.5,6.9+3.8-d,0,m2)

        #####################################################################################
        # GENERATE THE LINES
        #####################################################################################
        l1 = gmsh.model.geo.addLine(p1,p2)
        l2 = gmsh.model.geo.add_ellipse_arc(p2,c1,p4,p3)
        l3 = gmsh.model.geo.add_ellipse_arc(p3,c1,p4,p4)
        l4 = gmsh.model.geo.add_ellipse_arc(p4,c1,p4,p5)
        l5 = gmsh.model.geo.add_ellipse_arc(p5,c1,p4,p6)
        l6 = gmsh.model.geo.addLine(p6,p7)  
        l7 = gmsh.model.geo.addLine(p7,p8)
        l8 = gmsh.model.geo.addLine(p8,p9)
        l9 = gmsh.model.geo.add_ellipse_arc(p9,c2,p11,p10)
        l10 = gmsh.model.geo.add_ellipse_arc(p10,c2,p11,p11)
        l11 = gmsh.model.geo.add_ellipse_arc(p11,c2,p11,p12)
        l12 = gmsh.model.geo.add_ellipse_arc(p12,c2,p11,p13)
        l13 = gmsh.model.geo.addLine(p13,p14)
        l14 = gmsh.model.geo.addLine(p14,p15)
        l15 = gmsh.model.geo.addLine(p15,p1)
        l16 = gmsh.model.geo.addLine(p4,p21)
        l17 = gmsh.model.geo.addLine(p3,p17)
        l18 = gmsh.model.geo.addLine(p5,p19)
        l19 = gmsh.model.geo.addLine(p17,p16)
        l20 = gmsh.model.geo.addLine(p16,p12)
        l21 = gmsh.model.geo.addLine(p21,p20)
        l22 = gmsh.model.geo.addLine(p20,p11)
        l23 = gmsh.model.geo.addLine(p19,p18)
        l24 = gmsh.model.geo.addLine(p18,p10)
        l25 = gmsh.model.geo.addLine(p16,p20)
        l26 = gmsh.model.geo.addLine(p17,p21)
        l27 = gmsh.model.geo.addLine(p21,p19)
        l28 = gmsh.model.geo.addLine(p20,p18)

        #####################################################################################
        # GENERATE THE SURFACES
        #####################################################################################
        cloop1 = gmsh.model.geo.addCurveLoop([l1,l2,l17,l19,l20,l12,l13,l14,l15])
        surface1 = gmsh.model.geo.addPlaneSurface([cloop1],1)
        cloop2 = gmsh.model.geo.addCurveLoop([l3,l16,-l26,-l17])
        surface2 = gmsh.model.geo.addPlaneSurface([cloop2],2)
        cloop3 = gmsh.model.geo.addCurveLoop([l4,l18,-l27,-l16])
        surface3 = gmsh.model.geo.addPlaneSurface([cloop3],3)
        cloop4 = gmsh.model.geo.addCurveLoop([l5,l6,l7,l8,l9,-l24,-l23,-l18])
        surface4 = gmsh.model.geo.addPlaneSurface([cloop4],4)
        cloop5 = gmsh.model.geo.addCurveLoop([l26,l21,-l25,-l19])
        surface5 = gmsh.model.geo.addPlaneSurface([cloop5],5)
        cloop6 = gmsh.model.geo.addCurveLoop([l27,l23,-l28,-l21])
        surface6 = gmsh.model.geo.addPlaneSurface([cloop6],6)
        cloop7 = gmsh.model.geo.addCurveLoop([l25,l22,l11,-l20])
        surface7 = gmsh.model.geo.addPlaneSurface([cloop7],7)
        cloop8 = gmsh.model.geo.addCurveLoop([l28,l24,l10,-l22])
        surface8 = gmsh.model.geo.addPlaneSurface([cloop8],8)

        #####################################################################################
        # VECTOR INCLUDING ALL THE LINES AND POINTS OF THE BOUNDARY
        #####################################################################################
        boundlines  = [l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15]
        boundpoints = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15]

        #####################################################################################
        # SYNCRONIZATION WITH THE KERNEL
        #####################################################################################
        gmsh.model.geo.synchronize()

        #####################################################################################
        # DEFINITION OF THE PHYSICAL GROUPS TO EXPORT NODES IN THE BOUNDARY. ONE GLOBAL, ONE PER LINE
        #####################################################################################
        gmsh.model.addPhysicalGroup(gdim-2, boundpoints,1) 
        gmsh.model.setPhysicalName(gdim-2, 1, "Boundary_Points")
        gmsh.model.addPhysicalGroup(gdim-1, boundlines,1) 
        gmsh.model.setPhysicalName(gdim-1, 1, "Boundary_Lines")
        for (i,boundline) in enumerate(boundlines):
            gmsh.model.addPhysicalGroup(gdim-1, [boundline],i+2) 
            gmsh.model.setPhysicalName(gdim-1, 1, "Boundary_L{}".format(i+1))

        #####################################################################################
        # DEFINITION OF THE 2D MESH
        #####################################################################################
        gmsh.model.mesh.generate(gdim)

        #####################################################################################
        # NODES IN THE BOUNDARY (THE GEOMETRICAL ONES AND THE CREATED ONES WHEN MESHING)
        #####################################################################################
        nodes_boundary = []
        for tag in gmsh.model.getEntitiesForPhysicalGroup(gdim-2,1):
            nodeTags, _, _ = gmsh.model.mesh.getNodes(gdim-2, tag)
            nodes_boundary.extend(nodeTags-1)
        for tag in gmsh.model.getEntitiesForPhysicalGroup(gdim-1,1):
            nodeTags, _,_ = gmsh.model.mesh.getNodes(gdim-1, tag)
            # In this case the label is "+1" because we need to jumpt two construction points created for the elliptical lines
            # (c1, c2) in the shape of the specimen. If in the geometry there are not construction poitns, put nodeTags-1
            nodes_boundary.extend(nodeTags+7) 
            
        #####################################################################################
        # SAVE THE POINTS IN A .DAT FILE
        #####################################################################################
        dat_filename = mesh_output+"/boundnodes/nodes_boundary.dat"
        with open(dat_filename, 'w') as dat_file:
            for node_tag in nodes_boundary:
                dat_file.write("{}\n".format(node_tag))

        #####################################################################################
        # OBTAIN THE NODES THAT ARE IN EACH PHYSICAL GROUP
        #####################################################################################
        numphys = 1 + len(boundlines) # number of physical groups
        
        # Nodes in the lines of the boundary (geometrical and created nodes when meshing)
        for i in range(2,numphys+1):
            # Create an empty vector to include the points 
            nodes_boundary = []
            # Include the initial point of the line
            nodeTags, _, _ = gmsh.model.mesh.getNodes(gdim-2, gmsh.model.getEntitiesForPhysicalGroup(gdim-2,1)[i-2])
            nodes_boundary.extend(nodeTags-1)
            # Include intermediate points
            for tag in gmsh.model.getEntitiesForPhysicalGroup(gdim-1,i):
                nodeTags, _,_ = gmsh.model.mesh.getNodes(gdim-1, tag)
                nodes_boundary.extend(nodeTags+7) # effect of construction points (look previous comments)
            # Include the final point of the line
            if i == numphys: # Include the initial point in the last line
                nodeTags, _, _= gmsh.model.mesh.getNodes(gdim-2, gmsh.model.getEntitiesForPhysicalGroup(gdim-2,1)[0])
            else:
                nodeTags, _, _= gmsh.model.mesh.getNodes(gdim-2, gmsh.model.getEntitiesForPhysicalGroup(gdim-2,1)[i-1])
            nodes_boundary.extend(nodeTags-1)
            # Save the nodes in a file .dat
            dat_filename = mesh_output+"/boundnodes/nodes_boundary_l{}.dat".format(i-1)
            with open(dat_filename, 'w') as dat_file:
                for node_tag in nodes_boundary:
                    dat_file.write("{}\n".format(node_tag))
        
        # Include geometrical nodes (starting from 0)
        geometrical_points = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
        for elemento in reversed(geometrical_points):
            nodes_boundary.insert(0, elemento)
        

        # Reproduce the mesh
        # gmsh.fltk.run()

        # Export the mesh (ply2 format)
        gmsh.write(mesh_output+"/mesh.ply2")

    #................. Close GMSH
    gmsh.finalize()
        
    