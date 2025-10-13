#####################################################################################
# THIS CODE APPLIES THE PHASE FIELD MODEL DEVELOPED IN GOLAHMAR ET AL. (2023) FOR
# FATIGUE IN THE FENICSX SOFTWARE. 
#####################################################################################
"""
@authors: Sara Jimenez Alfaro & Emilio Martinez Paneda
"""
# Last update: 08/05/2025
import dolfinx
from mpi4py import MPI
import petsc4py 
from petsc4py import PETSc
import ufl
from SNES_solver import SNES_problem
import numpy as np
import dolfinx.fem.petsc
from math import *

def phasefield_solver(name_input,mesh,cells,facets, Parameters, coords_point):

    # Parallel computation: multiprocessors
    comm = MPI.COMM_WORLD

    #############################################################
    ## INITIAL PARAMETERS
    #############################################################
    E, nu, Gc, l0  = Parameters.get('E'), Parameters.get('nu'), Parameters.get('Gc'), Parameters.get('l0')
    a, b    = Parameters.get("a"), Parameters.get("b")
    qmax    = Parameters.get('Load')*3.8/6.9
    Ncycles = Parameters.get('NCycles')
    Nitermax= 100
    Tol     = 1.0e-7
    fdim    = mesh.geometry.dim
    kres    = 1.e-6
    lmbda3d = E*nu / ((1+nu)*(1-2*nu))
    mu      = E / 2 / (1+nu)
    lmbda   = 2*mu*lmbda3d/(lmbda3d+2*mu)  # Plane stress
    #############################################################
    ## DAMAGE MODEL
    #############################################################

    if Parameters.get('model') == 'AT2':
        cw = dolfinx.fem.Constant(mesh, 1/2.0)
        wp = lambda d: 2.0*d
        gp = lambda d: (1-d)**2 + kres
        
    else:
        cw = dolfinx.fem.Constant(mesh, 2/3.0)
        wp = lambda d: 1.0
        gp = lambda d: (1-d)**2 + kres
        Sc = sqrt(3/8*E*Gc/l0)

    #############################################################
    ## ADDITIONAL FATIGUE PARAMETERS
    #############################################################
    qe  = a*1e6**(-b)*3.8/6.9
    n   = 0.5*1/b - 0.13
    Npe = 6.e5
    Spe = a*Npe**(-b)
    a0  = Npe * (Spe/Sc)**(2*n) / (1 - Spe/Sc)
    Ratio, kappa = -1, 0.5
    a_n = 1/2*Sc**2/E
    
    #############################################################
    ## FUNCTION SPACES
    #############################################################
    
    # Element for the history variable (calculation)
    Element_Scalar_hq = ufl.FiniteElement("Quadrature", mesh.ufl_cell(), degree=1, quad_scheme='default')

    # function space for the displacement
    V_u = dolfinx.fem.VectorFunctionSpace(mesh, ("CG",1))
    # function space for the damage 
    V_d = dolfinx.fem.FunctionSpace(mesh, ("CG",1))
    # function space for the history degradation
    V_h = dolfinx.fem.FunctionSpace(mesh, Element_Scalar_hq)
    # function space for the fatigue degradation
    V_f = dolfinx.fem.FunctionSpace(mesh, ("CG",1))

    # Functions
    # Functions for displacement and damage
    u = dolfinx.fem.Function(V_u, name='Displacement')
    d = dolfinx.fem.Function(V_d, name='Damage')
    f = dolfinx.fem.Function(V_f, name='Fatigue_Degradation')
    a = dolfinx.fem.Function(V_f, name='Fatigue_dF')
    a_e = dolfinx.fem.Function(V_f, name='Endurance_Energy')  
    # Old functions
    H_old = dolfinx.fem.Function(V_h, name='Old_EHistory')
    d_old = dolfinx.fem.Function(V_d, name='Old_Damage')  
    a_old = dolfinx.fem.Function(V_f, name='Old_FdF')   
    amax_old = dolfinx.fem.Function(V_f, name='Old_Max_FdF')  
    #............... Test functions
    utest = ufl.TestFunction(V_u)
    dtest = ufl.TestFunction(V_d)
    # Definition of the differential for integration
    ds  = ufl.Measure("ds", domain=mesh, subdomain_data=facets)
    dx  = ufl.Measure("dx", domain=mesh)

    #############################################################
    ## BOUNDARY CONDITIONS
    #############################################################

    # Definition of the load vector
    cycles = np.linspace(start=1, stop=Ncycles, num=Ncycles)
    Uo1 = dolfinx.fem.Constant(mesh, 0.0)
    Uo2 = dolfinx.fem.Constant(mesh, petsc4py.PETSc.ScalarType((0.0,0.0)))
    
    left_dofs1 = dolfinx.fem.locate_dofs_topological(V_u.sub(0), facets.dim, facets.indices[facets.values == 3])
    def point(x):
        return np.logical_and(np.isclose(x[0], coords_point[0]), np.isclose(x[1], coords_point[1]))
    left_dofs2 = dolfinx.fem.locate_dofs_topological(V_u, facets.dim-1, dolfinx.mesh.locate_entities_boundary(mesh,0,point))
    # .... Dirichlet boundary conditions displacement
    bc_left1  = dolfinx.fem.dirichletbc(Uo1, left_dofs1, V_u.sub(0))
    bc_left2  = dolfinx.fem.dirichletbc(Uo2, left_dofs2, V_u)
    # .... Storage of the dirichlet boundary conditions
    bcs       = [bc_left1,bc_left2]

    # Dirichlet bcs for damage
    right_dofs = dolfinx.fem.locate_dofs_topological(V_d, facets.dim, facets.indices[facets.values == 4])
    left_dofs  = dolfinx.fem.locate_dofs_topological(V_d, facets.dim, facets.indices[facets.values == 3])
    bc_right_d = dolfinx.fem.dirichletbc(0.0, right_dofs, V_d)
    bc_left_d  = dolfinx.fem.dirichletbc(0.0, left_dofs, V_d)
    bcd = [bc_right_d,bc_left_d]

    # Neumann boundary conditions and volume forces for displacements (2d)
    q = dolfinx.fem.Constant( mesh, petsc4py.PETSc.ScalarType((0.0,0.0)) )
    b = dolfinx.fem.Constant( mesh, petsc4py.PETSc.ScalarType((0.0,0.0)) )

    #############################################################
    ## COMPLEMENTARY FUNCTIONS IN THE PROBLEM
    #############################################################

    # Strain tensor
    def eps(v):
        return ufl.sym(ufl.grad(v))
        
    #... Non-damaged positive strain energy density
    def Psip_0(v):
        trace, det  = ufl.tr(eps(v)), ufl.det(eps(v))
        eps3 = -nu/(1-nu) * trace # tercera componente del tensor de deformaciones
        e1   = ufl.max_value(0.5 * ( trace + ufl.sqrt(trace**2-4*det)),eps3)
        e2   = ufl.max_value(0.5 * ( trace - ufl.sqrt(trace**2-4*det)),ufl.min_value(0.5 * ( trace + ufl.sqrt(trace**2-4*det)),eps3))
        e3   = ufl.min_value(0.5 * ( trace - ufl.sqrt(trace**2-4*det)),ufl.min_value(0.5 * ( trace + ufl.sqrt(trace**2-4*det)),eps3))
        psip0 = ufl.conditional(e3>0,0.5 * lmbda3d * (e1 + e2 + e3)**2 + mu * (e1**2 + e2**2 + e3**2),\
        ufl.conditional(e2 + nu*e3>0,0.5 * lmbda3d * (e1 + e2 + 2*nu*e3)**2 + mu * ((e1+nu*e3)**2 + (e2+nu*e3)**2),\
        ufl.conditional((1-nu)*e1+nu*(e2+e3)>0,0.5 * lmbda3d *(1+nu)/(nu*(1-nu**2)) * ((1-nu)*e1+nu*e2+nu*e3)**2,0)))
        return psip0
    #... Stress tensor
    def sigma(v,p):
        trace  = ufl.tr(eps(v))
        return  gp(p) * (lmbda * trace * ufl.Identity(fdim) + 2*mu * eps(v)) 
    # History variable
    if Parameters.get('model') == 'AT2':
        def H_function(v,Hold):
            return ufl.max_value(Psip_0(v),Hold)
    else:
        def H_function(v,Hold):
            return ufl.max_value(ufl.max_value(Psip_0(v),Hold),3*Gc/16/l0)
        

    #############################################################
    ## DISPLACEMENT PROBLEM
    #############################################################
    # Definition of the residual
    u_residual = ufl.inner(sigma(u,d_old),eps(utest)) * ufl.dx - ufl.inner(b, utest) * ufl.dx - ufl.inner(q, utest)*ds(4)
    # SNES solver for displacement
    u_problem  = SNES_problem(u_residual, u, bcs)
    # Matrix and vector for displacement system
    F_u = dolfinx.la.create_petsc_vector(V_u.dofmap.index_map, V_u.dofmap.index_map_bs)
    J_u = dolfinx.fem.petsc.create_matrix(u_problem.J)
    # Solver for displacement problem
    u_solver = PETSc.SNES().create()
    u_solver.setType("ksponly") 
    u_solver.setTolerances(rtol=1.0e-10, max_it=50)
    # Linear system solver
    u_solver.getKSP().setType("preonly")
    u_solver.getKSP().setTolerances(rtol=1.0e-10)
    u_solver.getKSP().getPC().setType("lu")
    # Variables in the Newton problem
    u_solver.setFunction(u_problem.Residual, F_u)
    u_solver.setJacobian(u_problem.Jacobian, J_u)
    u_solver.getKSP().getPC().setFactorSolverType('mumps')


    # Fatigue driving force function (f2)
    Heaviside = lambda x: ufl.conditional(x<0, 0.0, 1.0)  
    def Heaviside(x):
        return ufl.conditional(x<0,0.0,1.0)
    def f_function(a):
        return ufl.conditional(a>=0, ufl.conditional(a<=a0, (1-a/a0)**2,0.0), 1.0)
    def amax_function(amaxT,v):
        return ufl.max_value(amaxT,Psip_0(v))
    def a_function(alpha,amaxT,v,a_e):
        amax = amax_function(amaxT,v)
        Delta_a = (amax/a_n)**n * ((1-Ratio)/2)**(2*kappa*n) * Heaviside(amax*((1-Ratio)/2)**(2*kappa) - a_e)
        return alpha + Delta_a 

    #############################################################
    ## DAMAGE PROBLEM
    #############################################################

    # Damage residual with the gradient term
    term1 = -2* (1-d) * H_function(u,H_old) * dtest
    term2 = Gc*l0/(2*cw) * f * ( wp(d)/(2*l0**2) * dtest + ufl.inner(ufl.grad(d),ufl.grad(dtest)) )
    term3 = - Gc*l0/(2*cw) * ufl.inner( ufl.grad(f),ufl.grad(d) ) * dtest
    d_residual = ( term1 + term2 + term3 ) * ufl.dx     # SNES solver for damage
    d_problem  = SNES_problem(d_residual, d, bcd)
    # Matrix and vector for damage system
    F_d = dolfinx.la.create_petsc_vector(V_d.dofmap.index_map, V_d.dofmap.index_map_bs)
    J_d = dolfinx.fem.petsc.create_matrix(d_problem.J)
    # Solver for damage problem
    d_solver = PETSc.SNES().create()
    d_solver.setType("vinewtonrsls") 
    d_solver.setTolerances(atol=1.0e-10)
    d_solver.setTolerances(rtol=1.0e-10, max_it=50)
    # We define the linear system solver
    d_solver.getKSP().setType("preonly")
    d_solver.getKSP().setTolerances(atol=1.0e-10)
    d_solver.getKSP().setTolerances(rtol=1.0e-10)
    d_solver.getKSP().getPC().setType("lu")
    # We define the variables of the Newton problem
    d_solver.setFunction(d_problem.Residual, F_d)
    d_solver.setJacobian(d_problem.Jacobian, J_d)
    d_solver.getKSP().getPC().setFactorSolverType('mumps')

    # Damage limits
    d_up  = dolfinx.fem.Function(V_d, name='Upper_Bound_Damage')
    d_low = dolfinx.fem.Function(V_d, name='Lower_Bound_Damage')
    # Include bcd
    dolfinx.fem.set_bc(d_up.vector,bcd)
    dolfinx.fem.set_bc(d_low.vector,bcd)
    # Define the bounds
    with d_up.vector.localForm() as bc_local:
        bc_local.set(1.0)
    with d_low.vector.localForm() as bc_local:
        bc_local.set(-0.4)
    # We set the bound (Note: they are passed as reference and not as values)
    d_solver.setVariableBounds(d_low.vector,d_up.vector)

    ################################################################
    # OUTCOME DOCUMENT
    ################################################################
    
    from pathlib import Path
    results_folder = Path(str("output/profiles/"))
    ffiled = dolfinx.io.XDMFFile(mesh.comm, results_folder / str(name_input+"_output.xdmf"), "w")
    ffiled.write_mesh(mesh)
    if comm.rank==0:
        Numerics = []

    
    ############################################################
    ## INITIAL PROBLEM FOR THE ENDURANCE ENERGY
    #############################################################

    q.value = petsc4py.PETSc.ScalarType((qe,0.0))
    # Displacement problem
    u_solver.solve(None, u.vector)
    u_solver.destroy, J_u.destroy(), F_u.destroy()
    u.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
    # Update the history variable
    a_e.interpolate(dolfinx.fem.Expression(Psip_0(u), V_f.element.interpolation_points()))
    a_e.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)  
    u.x.array[:] = 0.
    ffiled.write_function(a_e)
    

    #................ Total area/volume in the solid
    vol_rank  = dolfinx.fem.assemble_scalar(dolfinx.fem.form(dolfinx.fem.Constant(mesh, 1.0) * ufl.dx))
    vol_total = MPI.COMM_WORLD.allreduce(vol_rank, op=MPI.SUM)  

    mesh_global = np.concatenate(comm.allgather(np.array(mesh.geometry.x)))
    mesh_global = np.array(list({tuple(row) for row in mesh_global}))
    midpoint_position = np.where((mesh_global[:, 0] == 32.5) & (mesh_global[:, 1] == 6.9))[0]
    
    #############################################################
    ## INTERMEDIATE INTERPOLATIONS
    #############################################################
    import basix
    class interpolate_quadrature:
        """
        class for quadrature level projections: suitable for quadrature elements
        """
        def __init__(self, ufl_expr, fem_func:dolfinx.fem.Function):
            q_dim = fem_func.function_space._ufl_element.degree()
            mesh = fem_func.ufl_function_space().mesh
            basix_celltype = basix.cell.string_to_type(mesh.ufl_cell().cellname())
            quadrature_points, weights = basix.make_quadrature(basix_celltype, q_dim)
            map_c = mesh.topology.index_map(mesh.topology.dim)
            num_cells = map_c.size_local + map_c.num_ghosts
            cells = np.arange(0, num_cells, dtype=np.int32)

            self.mesh = mesh
            self.cells =  cells
            self.quadrature_points =  quadrature_points
            self.ufl_expr = ufl_expr
            self.fem_func = fem_func

        def eval(self):
            expr_expr = dolfinx.fem.Expression(self.ufl_expr, self.quadrature_points)
            expr_eval = expr_expr.eval(self.mesh, self.cells)
            self.fem_func.x.array[:] = expr_eval.flatten()[:] 
    
    # Mapping of the history variable
    map_Hold = interpolate_quadrature(H_function(u,H_old),H_old)
    H_old.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
    
    #############################################################
    ## ITERATION
    #############################################################
    
    import time
    t0 = time.time()
    flag_break = 0

    for nc in cycles[0:]:

        niter = 0

        q.value = petsc4py.PETSc.ScalarType((qmax,0.0))

        for niter in range(Nitermax):

            # Displacement problem
            u_solver.solve(None, u.vector)
            u_solver.destroy, J_u.destroy(), F_u.destroy()
            u.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

            # Update the degradation function
            a.interpolate(dolfinx.fem.Expression(a_function(a_old,amax_old,u,a_e), a.function_space.element.interpolation_points()))
            a.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
            f.interpolate(dolfinx.fem.Expression(f_function(a), f.function_space.element.interpolation_points()))
            f.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

            # Damage problem
            d_solver.solve(None, d.vector)
            d_solver.destroy, J_d.destroy(), F_d.destroy()
            d.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

            if d_solver.getConvergedReason()<0 and comm.rank == 0:
                print('Damage solver')
                print("Number of iterations : ", d_solver.getIterationNumber())
                print(f"Converged reason = {d_solver.getConvergedReason():1.3f}")
                print(f"Function norm = {d_solver.getFunctionNorm():.4e}")
                # ffiled.close()
                # exit()
         
            # Check convergence
            error_rank  = dolfinx.fem.assemble_scalar(dolfinx.fem.form( ufl.inner( d-d_old,d-d_old )*ufl.dx ))
            error_total = MPI.COMM_WORLD.allreduce(error_rank, op=MPI.SUM)
            errorL2     = np.sqrt(error_total)/vol_total

            # Update damage for the next iteration
            d.vector.copy(d_old.vector)
            d_old.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

            # print(errorL2)

            if errorL2 < Tol:
                break
        
        # Update the history variable
        map_Hold.eval()
        H_old.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)  

        # Update the fatigue driving force
        a.vector.copy(a_old.vector)
        a_old.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

        # Update the energy term in the history
        amax_old.interpolate(dolfinx.fem.Expression(amax_function(amax_old,u), V_f.element.interpolation_points()))
        amax_old.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

        if comm.rank==0 and nc%Parameters.get('Njump')==0:
            Numerics.append([nc,niter,errorL2])
            np.savetxt(str("output/" + str( name_input + '_numerics.dat')), Numerics)

        if nc%Parameters.get('Njump')==0:
            # Save in a xdmf file
            ffiled.write_function(d,nc)
            ffiled.write_function(u,nc)
            ffiled.write_function(a,nc)
            ffiled.write_function(f,nc)

        d_global = np.concatenate(comm.allgather(np.array(d.x.array)))
        
        if d_global[midpoint_position] >= 0.95:
            
            # This loop captures when the specimen is broken (when at the center the phase field length scale is bigger than 0.95)
            flag_break = 1
            if comm.rank==0:
                Numerics.append([nc,niter,errorL2])
                np.savetxt(str("output/" + str( name_input + '_numerics.dat')), Numerics)
                
        if flag_break == 1:
            break

    ffiled.close()
    tf = time.time()
    if comm.rank == 0:
        print('time (min)',(tf-t0)/60)
                    