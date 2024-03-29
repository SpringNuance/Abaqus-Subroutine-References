# Phase field fracture implementation in FEniCS    
# The code is distributed under a BSD license     
      
# If using this code for research or industrial purposes, please cite:
# Hirshikesh, S. Natarajan, R. K. Annabattula, E. Martinez-Paneda.
# Phase field modelling of crack propagation in functionally graded materials.
# Composites Part B: Engineering 169, pp. 239-248 (2019)
# doi: 10.1016/j.compositesb.2019.04.003
      
# Emilio Martinez-Paneda (mail@empaneda.com)
# University of Cambridge

# Preliminaries and mesh
from dolfin import *
mesh = Mesh('mesh.xml')

# Define Space
V = FunctionSpace(mesh, 'CG', 1)
W = VectorFunctionSpace(mesh, 'CG', 1)
WW = FunctionSpace(mesh, 'DG', 0)
p, q = TrialFunction(V), TestFunction(V)
u, v = TrialFunction(W), TestFunction(W)

# Introduce manually the material parameters of the two compounds
K = 0.2;
E2 = 380e3;
E1 = 210e3;
nu2 = 0.26;
nu1 = 0.31;
KIc2 = 5.2;  #MPasqrt(m)
KIc1 = 9.6;


class Mori_Tanaka:
    def __init__(self, material_parameters,mesh):
        mp = material_parameters
        self.mesh = mesh
        self.K = mp['K']
        self.E2 = mp['E2']
        self.E1 = mp['E1']
        self.nu2 = mp['nu2']
        self.nu1 = mp['nu1']
        self.KIc2 = mp['KIc2']
        self.KIc1 = mp['KIc1']

    def getFunctionMaterials(self, V):
        self.x = SpatialCoordinate(self.mesh)
        self.Vf = pow((0.5 + self.x[1]), self.K) 
        self.K1 = self.E1/(3.0*(1.-2.0*self.nu1))
        self.K2 = self.E2/(3.0*(1.-2.0*self.nu2))
        self.G1 = self.E1/(2.0*(1.0+ self.nu1))
        self.G2 = self.E2/(2.0*(1.0+ self.nu2))
        self.Ke = self.K1 + (self.K2-self.K1)*(1.-self.Vf)\
		/(1.+3.*self.Vf*((self.K2-self.K1)/(3.*self.K1 + 4.*self.G1)))
        term1 = self.Vf*(self.G2 - self.G1);
        term2 = (9.0*self.K1+8.0*self.G1)/(6.0*(self.K1+2.0*self.G1))
        self.Ge = self.G1 + (self.G2-self.G1)*(1.-self.Vf)\
		/(1.0+term1/(self.G1+self.G1*term2))      
        self.E = 9.0*self.Ke*self.Ge/(3.0*self.Ke + self.Ge)
        self.nu = (3.0*self.Ke - 2.0*self.Ge)/(2.0*( 3.0*self.Ke + self.Ge ))
        term3 = self.E/(1.0 - self.nu**2)
        term4 = self.Vf*( ((1.-self.nu1**2)/self.E1)*(self.KIc1/self.KIc2)**2 )
        term5 = (1.- self.Vf)*((1.-self.nu2**2)/self.E2)
        self.KIc = self.KIc2*sqrt(  term3*(term4 + term5)  )
        fac = (1.- self.nu**2 )*1e3 # for plane-strain
        self.Gc = fac*( self.KIc**2/self.E )
        effectiveMdata = {'E':self.E, 'nu':self.nu, 'Gc':self.Gc}
        return effectiveMdata

material_parameters = {'K':K, 'E2':E2, 'E1':E1, 'nu2':nu2, 'nu1': nu1, 'KIc2': KIc2, 'KIc1':KIc1}
mat = Mori_Tanaka(material_parameters, mesh)
EmatData = mat.getFunctionMaterials(V)

E  = EmatData['E']
nu = EmatData['nu']
Gc = EmatData['Gc']

lmbda, mu = (E*nu/((1.0 + nu )*(1.0-2.0*nu))) , (E/(2*(1+nu)))
l = 0.015


# Constituive functions
def epsilon(u):
    return sym(grad(u))
def sigma(u):
    return 2.0*mu*epsilon(u)+lmbda*tr(epsilon(u))*Identity(len(u))
def psi(u):
    return 0.5*(lmbda+mu)*(0.5*(tr(epsilon(u))+abs(tr(epsilon(u)))))**2+\
           mu*inner(dev(epsilon(u)),dev(epsilon(u)))		
def H(uold,unew,Hold):
    return conditional(lt(psi(uold),psi(unew)),psi(unew),Hold)
		
# Boundary conditions
top = CompiledSubDomain("near(x[1], 0.5) && on_boundary")
bot = CompiledSubDomain("near(x[1], -0.5) && on_boundary")
def Crack(x):
    return abs(x[1]) < 1e-03 and x[0] <= 0.0
load = Expression("t", t = 0.0, degree=1)
bcbot= DirichletBC(W, Constant((0.0,0.0)), bot)
bctop = DirichletBC(W.sub(1), load, top)
bc_u = [bcbot, bctop]
bc_phi = [DirichletBC(V, Constant(1.0), Crack)]
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundaries.set_all(0)
top.mark(boundaries,1)
ds = Measure("ds")(subdomain_data=boundaries)
n = FacetNormal(mesh)

# Variational form
unew, uold = Function(W), Function(W)
pnew, pold, Hold = Function(V), Function(V), Function(V)
E_du = ((1.0-pold)**2)*inner(grad(v),sigma(u))*dx
E_phi = (Gc*l*inner(grad(p),grad(q))+((Gc/l)+2.0*H(uold,unew,Hold))\
            *inner(p,q)-2.0*H(uold,unew,Hold)*q)*dx
p_disp = LinearVariationalProblem(lhs(E_du), rhs(E_du), unew, bc_u)
p_phi = LinearVariationalProblem(lhs(E_phi), rhs(E_phi), pnew, bc_phi)
solver_disp = LinearVariationalSolver(p_disp)
solver_phi = LinearVariationalSolver(p_phi)

# Initialization of the iterative procedure and output requests
t = 0
u_r = 0.007
deltaT  = 0.1
tol = 1e-3
conc_f = File ("./ResultsDir/phi.pvd")
fname = open('ForcevsDisp.txt', 'w')

# Staggered scheme
while t<=1.0:
    t += deltaT
    if t >=0.25:
        deltaT = 0.0001
    load.t=t*u_r
    iter = 0
    err = 1

    while err > tol:
        iter += 1
        solver_disp.solve()
        solver_phi.solve()
        err_u = errornorm(unew,uold,norm_type = 'l2',mesh = None)
        err_phi = errornorm(pnew,pold,norm_type = 'l2',mesh = None)
        err = max(err_u,err_phi)
        
        uold.assign(unew)
        pold.assign(pnew)
        Hold.assign(project(psi(unew), WW))

        if err < tol:
		
            print ('Iterations:', iter, ', Total time', t)

            if round(t*1e4) % 10 == 0:
                conc_f << pnew

                Traction = dot(sigma(unew),n)
                fy = Traction[1]*ds(1)
                fname.write(str(t*u_r) + "\t")
                fname.write(str(assemble(fy)) + "\n")
	    	    
fname.close()
print ('Simulation completed') 