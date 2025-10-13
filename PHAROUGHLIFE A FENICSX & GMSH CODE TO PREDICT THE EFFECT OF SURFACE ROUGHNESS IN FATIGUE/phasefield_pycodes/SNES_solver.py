#####################################################################################
# THIS CODE IS THE SNES SOLVER USED TO SOLVE THE PHASE FIELD MODEL
#####################################################################################
"""
@authors: Sara Jimenez Alfaro & Emilio Martinez Paneda
"""
# Last update: 08/05/2025
import dolfinx
from petsc4py import PETSc
import ufl

class SNES_problem():
    
    def __init__(self,F,u,bcs):
        
        """ Definition of the class atributes """
        du      = ufl.TrialFunction(u.function_space)      # du = increment solution
        self.F  = dolfinx.fem.form(F)                      # F = residual       
        self.bc = bcs                                      # bcs = boundary conditions
        self.u  = u                                        # u = solution 
        self.J  = dolfinx.fem.form(ufl.derivative(F,u,du)) # J = jacobian

        # .. Generation of matrices in the no-linear problem
        self.Jmatrix = dolfinx.fem.petsc.create_matrix(self.J)
        self.Fvector = dolfinx.fem.petsc.create_vector(self.F)
        
    def Residual(self, snes, x, Fvector):

        """Assemble residual vector."""

        # .... Connection between ranks in the vector x: FORWARD
        x.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
        # ..... We copy the vector x in the vector u
        x.copy(self.u.vector)
        # .... Connection between ranks in the vector u in the class
        self.u.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

        # Preallocate values 0.0 in the vector F (en todos los fantasma y los owner)
        with Fvector.localForm() as f_local:
            f_local.set(0.0)

        # Create the vector in a petsc way
        dolfinx.fem.petsc.assemble_vector(Fvector, self.F)
        # .... Include the boundary conditions using the lifting technique
        dolfinx.fem.petsc.apply_lifting(Fvector, [self.J], bcs=[self.bc], x0=[x], scale=-1.0)
         # .... Connection between ranks in the vector x: REVERSE
        Fvector.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        # .... Redefine the function F with the boundary condition.
        dolfinx.fem.petsc.set_bc(Fvector, self.bc, x, -1.0)

        Fvector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
        
    def Jacobian(self, snes, x, Jmatrix, P):
        
        """ Assemble the Jacobian matrix """
        # .......... Zero preallocation in the matrix
        Jmatrix.zeroEntries()
        #.......... Assemblage of the Jacobian. Remove columns and rows related to DirichletBC
        dolfinx.fem.petsc.assemble_matrix(Jmatrix, self.J, self.bc)
        Jmatrix.assemble()
        