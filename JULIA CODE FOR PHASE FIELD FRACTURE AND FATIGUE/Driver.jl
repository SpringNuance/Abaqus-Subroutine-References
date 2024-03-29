#Julia Code for Phase Field Fatigue
#If using this code please cite:
#P.K. Kristensen, A. Golahmar, E. Martinez-Paneda, C.F. Niordson. 
#Accelerated high-cycle phase field fatigue predictions. 
#European Journal of Mechanics A/Solids 100, 104991 (2023)

include("Subroutines.jl")
#Define problem domain. Here a 2D unit  square with 100x100 quad elements
grid = generate_grid(Quadrilateral, (100,100), Vec((-0.5, -0.5)), Vec((0.5, 0.5)));
#meshes can be read from gmsh or abaqus files as well, see the Ferrite.jl documentation for more information
#grid = saved_file_to_grid("SENT.msh")
#Define dimension of problem (2D)
dim=2;

#Parameters of timesteps and solver
n_timesteps = 1000;
u_max = 0.0008;
plot_frequency = 200;
history_output_frequency=100;

loadingtype = FatigueCLA #Monotonic, Fatigue, FatigueCLA

#Parameters for ModifiedNewtonRaphson
n_c  = 10;
nᵢ = 20;
rebuild_frequency_base = 3;


#\mathcal{N} in the paper, (FatigueCLA only)
CyclesPerIncrement=1.0

#Crack endpoint and direction (for output purposes)
a₀=0.0
CrackDir = 1
outputset= "top" #name of the set at which loading is applied (see below) for output purposes


ElementShape = RefCube #Select triangular of square elements (works for both 2D and 3D) RefCube or RefTetrahedron
ElementOrder = 1 #Linear elements
QuadratureOrder = 2 #second order quadrature
#Define interpolation
ip = Lagrange{dim,ElementShape,ElementOrder}()
#Define quadrature rule
qr = QuadratureRule{dim,ElementShape}(QuadratureOrder)
#Add machinery for cell (element) values
cellvalues_u=CellVectorValues(qr,ip);
cellvalues_ϕ=CellScalarValues(qr,ip);




#Set up constraint handlers - first define or assign node and face sets
addnodeset!(grid,"Crack",x -> x[1]<=0.0 && (x[2]≈ 0 || x[2] ≈ 0.01)) #You can make any function of coordinates x to add the nodes you need to the crack set
addnodeset!(grid,"top",x -> x[2]≈ 0.5 ) 
∂Ω₁ = getnodeset(grid,"top") 
∂Ω₂ = getfaceset(grid,"bottom")#These are default sets in the mesh formed on line 4
∂Ω₃ = getnodeset(grid,"Crack")

#Define boundary conditions
dbc₁=Dirichlet(:u, ∂Ω₁, (x,t) -> t, 2) #Faceset "top" is prescribed displacement t on displacement dof 2
dbc₂ = Dirichlet(:u,∂Ω₂, (x,t) -> [0,0], [1,2]) #Faceset "bottom" is prescribed displacement 0 on displacement dofs 1 and 2
dbc₃ = Dirichlet(:ϕ,∂Ω₃, (x,t) -> 1) # nodeset "Crack" is prescribed the value 1 on the phase field dof

##################### THIS SECTION DOES NOT REQUIRE ALTERATION ################
#Set up dofHandler                                                            #
dh_u=DofHandler(grid)                                                         #
#Add a dof named u with dimension 2 (vector)                                  #
push!(dh_u,:u,dim)                                                            #
close!(dh_u);                                                                 #
dh_phi=DofHandler(grid)                                                       #
#Add a dof named ϕ with dimension 1 (scalar)                                  #
push!(dh_phi,:ϕ,1)                                                            #
close!(dh_phi);                                                               #
#Define constraint handlers and add boundary conditions to them               #
ch_u = ConstraintHandler(dh_u);                                               #
add!(ch_u,dbc₁)                                                               #
add!(ch_u,dbc₂)                                                               #
close!(ch_u)                                                                  #
update!(ch_u,0.0);                                                            #
ch_phi = ConstraintHandler(dh_phi);                                           #
add!(ch_phi,dbc₃)                                                             #
close!(ch_phi)                                                                #
update!(ch_phi,0.0);                                                          #
###############################################################################

#Define material properties
Material = Brittle(210000.,0.3,2.7,0.04,Isotropic,loadingtype,CyclesPerIncrement, dim);
#E, ν, Gc, ℓ, decomposition flag, Loading type, Cycles pr increment
#Decomposition flag options: Isoptropic, VolDev, Spectral, NoTension

#Set up output struct
output=OutputVariables(plot_frequency,history_output_frequency,a₀,CrackDir,outputset);

#Set up solver state struct
solver=SolverState(rebuild_frequency_base,n_c,u_max,n_timesteps, nᵢ,2000,1.e-5,1.e-4,loadingtype);

#Start solving problem
@time Problem(solver,output, cellvalues_ϕ,cellvalues_u, dh_phi, dh_u,
           ch_phi,ch_u, grid, Material);

