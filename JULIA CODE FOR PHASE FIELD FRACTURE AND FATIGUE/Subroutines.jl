#Julia Code for Phase Field Fatigue
#If using this code please cite:
#P.K. Kristensen, A. Golahmar, E. Martinez-Paneda, C.F. Niordson. 
#Accelerated high-cycle phase field fatigue predictions. 
#European Journal of Mechanics A/Solids 100, 104991 (2023)

#Load necessary packages
using Ferrite, SparseArrays, Tensors, Printf, LinearAlgebra

#define some useful enumerators
@enum SolutionType NewtonRaphson ModifiedNewtonRaphson
@enum Load Monotonic Fatigue FatigueCLA
@enum StrainDecomp Isotropic VolDev Spectral NoTension

#Data strucutre for material
struct Brittle{T, S <: SymmetricTensor{4, 3, T}}
    G::T  # Shear modulus
    K::T  # Bulk modulus
    Gc::T # Fracture Toughness
    ℓ::T  # Phase field length scale
    flag::StrainDecomp
    Dᵉ::S # Elastic stiffness tensor
    load::Load
    Ncyc::T
    dim::Int64
end;

#Data strucutre for history variables
struct MaterialState{T}
    # Store history values
    H::T # History variable
    ϕ::T #phase field variable from last increment
    α::T #AccumulatedFatigue
    ψ::T #strain energy from last increment (fatigue only)
end

#Function for filling out the material struct
function Brittle(E, ν,Gc, ℓ, flag, load, Ncyc, dim)
    δ(i,j) = i == j ? 1.0 : 0.0 # helper function
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)
    temp(i,j,k,l) = 2.0G *( 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) + ν/(1.0-2.0ν)*δ(i,j)*δ(k,l))
    Dᵉ = SymmetricTensor{4, 3}(temp)
    return Brittle(G, K, Gc, ℓ, flag, Dᵉ,load, Ncyc, dim)
end;

#Initialize empty material state
function MaterialState()
    return MaterialState(0.0, 0.0, 0.0, 0.0)
end

#Data structure for solver information
mutable struct SolverState{T,F}
    rebuild_counter::T
    rebuild_frequency_base::T
    rebuild_frequency_min::T
    rebuild_start::Bool
    no_restart::T
    strategy::SolutionType
    u_max::F
    timesteps::T
    nitr_inner::T
    nitr_outer::T
    NewtonTOL_inner::F
    NewtonTOL_outer::F
    loading::Load
end

#Initialize data structure for solver information
function SolverState(RebuildFrequencyBase,RebuildPrN,Umax,nT,nitr_inner,nitr_outer,TOL_inner,TOL_outer,load)
    return SolverState(0,RebuildFrequencyBase,RebuildPrN,true,0,NewtonRaphson,Umax,nT,nitr_inner,nitr_outer,TOL_inner,TOL_outer,load)
end

#Data structure for output
mutable struct OutputVariables{T}
    plotframe::T
    totalIterations_outer::T
    totalIterations_phi::T
    totalIterations_u::T
    matrixFactorizations_u::T
    matrixFactorizations_phi::T
    plotFrequency::T
    historyFrequency::T
    a0::Float64
    CrackDir::T
    OutputSet::String
end
#Initialize data structure for output
function OutputVariables(field_frequency,history_frequency,a0,CrackDir,outputset)
    return OutputVariables(0,0,0,0,0,0,field_frequency,history_frequency,a0,CrackDir,outputset)
end
#Element stiffness matrix - displacement
function assemble_element_u!(Ke::Matrix, fe::Vector, cellvalues_u::CellVectorValues,
        ue::Vector,material::Brittle,state,state_old)
    nbase_u = getnbasefunctions(cellvalues_u)
    #Loop over integration points
    D = material.Dᵉ
    for q_point in 1:getnquadpoints(cellvalues_u)
        #dvol,
        dΩᵤ=getdetJdV(cellvalues_u,q_point)
        #Total strain
        if material.dim == 2
            ε_PlaneStrain=function_symmetric_gradient(cellvalues_u,q_point,ue)
            ∇ˢu = SymmetricTensor{2,3,Float64}((ε_PlaneStrain[1,1],ε_PlaneStrain[1,2], 0.0,
                ε_PlaneStrain[2,2], 0.0, 0.0))
        elseif material.dim == 3
            ∇ˢu =function_symmetric_gradient(cellvalues_u,q_point,ue)
        else
            error("Invalid element dimension")
        end
        #Phase field value
        ϕ = state[q_point].ϕ
        #(undegraded) Stress
        σ=D ⊡ ∇ˢu
        #Strain energy
        Ψ = CrackDrive(σ, ∇ˢu, material)
        #Recover state vairables
        ϕₙ = state_old[q_point].ϕ
        Hₙ = state_old[q_point].H
        αₙ = state_old[q_point].α
        ψₙ = state_old[q_point].ψ
        #Ensure phase field irreversibillity
        if material.load == Monotonic
            Δα = 0.
        elseif material.load ==Fatigue
            Δα = Ψ > ψₙ ? (Ψ-Ψₙ) : 0
        elseif material.load == FatigueCLA
            Δα = Ψ*material.Ncyc
        end
        α = αₙ + Δα
        #Update state variables
        H = max(Ψ,Hₙ)
        state[q_point]=MaterialState(H,ϕ,α,Ψ)
        #Phase field degradation based on old phi
        gdn= (1.0-ϕ)^2 +1.e-7

        ##ASSEMBLE DISPLACEMENT PROBLEM
        #Loop over displacement test functions
        for i in 1:nbase_u

            if material.dim == 2
                δε_2d = shape_symmetric_gradient(cellvalues_u,q_point,i)
                δε = SymmetricTensor{2,3,Float64}((δε_2d[1,1],δε_2d[1,2], 0.0,
                    δε_2d[2,2], 0.0, 0.0))
            else
                δε =shape_symmetric_gradient(cellvalues_u,q_point,i)
            end
            #Add contribution to element rhs

            fe[i] += gdn*(δε ⊡ σ)*dΩᵤ
            #Loop over displacement  trial shape functions
            for j in 1:i
                if material.dim ==2
                    ε_2d = shape_symmetric_gradient(cellvalues_u,q_point,j)
                    ε = SymmetricTensor{2,3,Float64}((ε_2d[1,1],ε_2d[1,2], 0.0,
                        ε_2d[2,2], 0.0, 0.0))
                else
                    ε = shape_symmetric_gradient(cellvalues_u,q_point,j)
                end
                Ke[i,j] += gdn*δε ⊡ D ⊡ ε * dΩᵤ
            end
        end
        ##ASSEMBLE PHASE FIELD PROBLEM
    end
    symmetrize_lower!(Ke)
    return Ke, fe
end
#Element residual - displacement
function assemble_residual_u!(fe::Vector, cellvalues_u::CellVectorValues,
        ue::Vector,material::Brittle,state,state_old, store::Bool)
    nbase_u = getnbasefunctions(cellvalues_u)
    #Loop over integration points
    D = material.Dᵉ
    for q_point in 1:getnquadpoints(cellvalues_u)
        #dvol,
        dΩᵤ=getdetJdV(cellvalues_u,q_point)
        #Total strain
        if material.dim == 2
            ε_PlaneStrain=function_symmetric_gradient(cellvalues_u,q_point,ue)
            ∇ˢu = SymmetricTensor{2,3,Float64}((ε_PlaneStrain[1,1],ε_PlaneStrain[1,2], 0.0,
                ε_PlaneStrain[2,2], 0.0, 0.0))
        elseif material.dim == 3
            ∇ˢu =function_symmetric_gradient(cellvalues_u,q_point,ue)
        else
            error("Invalid element dimension")
        end
        #Phase field value
        ϕ = state[q_point].ϕ
        #(undegraded) Stress
        σ=D ⊡ ∇ˢu
        #Strain energy
        Ψ = CrackDrive(σ, ∇ˢu, material)
          #Recover state vairables
          ϕₙ = state_old[q_point].ϕ
          Hₙ = state_old[q_point].H
          αₙ = state_old[q_point].α
          ψₙ = state_old[q_point].ψ
          #Ensure phase field irreversibillity
          if material.load == Monotonic
              Δα = 0.
          elseif material.load ==Fatigue
              Δα = Ψ > ψₙ ? (Ψ-Ψₙ) : 0
          elseif material.load == FatigueCLA
              Δα = Ψ*material.Ncyc
          end
          α = αₙ + Δα
          #Update state variables
          H = max(Ψ,Hₙ)
          if store
            state[q_point]=MaterialState(H,ϕ,α,Ψ)
          end


        #Phase field degradation based on old phi
        gdn= (1.0-ϕ)^2 +1.e-7

        ##ASSEMBLE DISPLACEMENT PROBLEM
        #Loop over displacement test functions
        for i in 1:nbase_u


            if material.dim == 2
                δε_2d = shape_symmetric_gradient(cellvalues_u,q_point,i)
                δε = SymmetricTensor{2,3,Float64}((δε_2d[1,1],δε_2d[1,2], 0.0,
                    δε_2d[2,2], 0.0, 0.0))
            else
                δε =shape_symmetric_gradient(cellvalues_u,q_point,i)
            end
            #Add contribution to element rhs

            fe[i] += gdn*(δε ⊡ σ)*dΩᵤ
        end
    end
    return fe
end
#Element stiffness matrix - phase field
function assemble_element_phi!(Ke::Matrix, fe::Vector, cellvalues_ϕ::CellScalarValues,
                                ϕe::Vector,material::Brittle,state,state_old)
    nbase_ϕ = getnbasefunctions(cellvalues_ϕ)

    Gc = material.Gc
    ℓ = material.ℓ

    # Loop over integration points
    for q_point in 1:getnquadpoints(cellvalues_ϕ)
        #dvol, should be the same, assuming they follow same interpolation
        dΩᵩ =getdetJdV(cellvalues_ϕ,q_point)
        #Phase field value
        ϕ = function_value(cellvalues_ϕ,q_point,ϕe)
        ∇ϕ = function_gradient(cellvalues_ϕ,q_point,ϕe)

        α=state[q_point].α
        H=state[q_point].H
        ψ=state[q_point].ψ
        #Recover state vairables
        Hₙ = state_old[q_point].H
        ϕₙ = state_old[q_point].ϕ
        ψₙ =state_old[q_point].ψ
        #Ensure phase field irreversibillity
        

        #Update state variables
        state[q_point]=MaterialState(H,ϕ,α,ψ)
        αₜ = Gc/12ℓ
        fdeg = α>αₜ ?   (2αₜ/(αₜ + α))^2 : 1.0
        #Phase field degradation function derivative

        gd′ =  -2.0(1.0-ϕ)
        ##ASSEMBLE PHASE FIELD PROBLEM
        for i in 1:nbase_ϕ
            δϕ = shape_value(cellvalues_ϕ,q_point,i)
            δ∇ϕ = shape_gradient(cellvalues_ϕ,q_point,i)
            fe[i] += (gd′*H*δϕ+fdeg*Gc/ℓ *δϕ * ϕ + fdeg*Gc*ℓ*δ∇ϕ ⋅ ∇ϕ)*dΩᵩ
            for j in 1:i
                ϕ′ = shape_value(cellvalues_ϕ,q_point,j)
                ∇ϕ′= shape_gradient(cellvalues_ϕ,q_point,j)
                gd′′ = 2.0ϕ′
                Ke[i,j] += (gd′′*H*δϕ+ fdeg*Gc/ℓ *δϕ*ϕ′+fdeg*Gc*ℓ * δ∇ϕ⋅∇ϕ′)*dΩᵩ
            end
        end
    end
    symmetrize_lower!(Ke)
    return Ke, fe
end
#Element residual - phase field
function assemble_residual_phi!(fe::Vector, cellvalues_ϕ::CellScalarValues,
                                ϕe::Vector,material::Brittle,state,state_old, store::Bool)
    nbase_ϕ = getnbasefunctions(cellvalues_ϕ)

    Gc = material.Gc
    ℓ = material.ℓ

    # Loop over integration points
    for q_point in 1:getnquadpoints(cellvalues_ϕ)
        #dvol, should be the same, assuming they follow same interpolation
        dΩᵩ =getdetJdV(cellvalues_ϕ,q_point)
        #Phase field value
        ϕ = function_value(cellvalues_ϕ,q_point,ϕe)
        ∇ϕ = function_gradient(cellvalues_ϕ,q_point,ϕe)

        α=state[q_point].α
        ψ=state[q_point].ψ
        #Recover state vairables
        Hₙ = state_old[q_point].H
        ϕₙ = state_old[q_point].ϕ
        ψₙ =state_old[q_point].ψ
        #Ensure phase field irreversibillity
        H=state[q_point].H

        #Update state variables
        if store
            state[q_point]=MaterialState(H,ϕ,α,ψ)
        end
        αₜ = Gc/12ℓ


        fdeg = α>αₜ ?   (2αₜ/(αₜ + α))^2 : 1.0


        #Phase field degradation function derivative
        gd′ =  -2.0(1.0-ϕ)

        ##ASSEMBLE PHASE FIELD PROBLEM
        for i in 1:nbase_ϕ
            δϕ = shape_value(cellvalues_ϕ,q_point,i)
            δ∇ϕ = shape_gradient(cellvalues_ϕ,q_point,i)
            fe[i] += (gd′*H*δϕ+fdeg*Gc/ℓ *δϕ * ϕ + fdeg*Gc*ℓ*δ∇ϕ ⋅ ∇ϕ)*dΩᵩ
        end
    end
    return fe
end
function symmetrize_lower!(K)
    for i in 1:size(K,1)
        for j in i+1:size(K,1)
            K[i,j] = K[j,i]
        end
    end
end;

#Subroutine for global assembly
function assemble_global(q::Vector, cellvalues,
                           K::SparseMatrixCSC,dh::DofHandler,
                           material::Brittle, states,states_old)
    #allocate element stiffness and rhs
    nbase = getnbasefunctions(cellvalues)

    Ke=zeros(nbase,nbase)
    fe = zeros(nbase)

    f=zeros(ndofs(dh))

    assembler = start_assemble(K,f)
    #THIS METHOD IS NOT ROBUST
    fielddim=nbase
    for (i, cell) in enumerate(CellIterator(dh))
        reinit!(cellvalues,cell)
        fill!(Ke,0)
        fill!(fe,0)
        eldofs=celldofs(cell)
        qe=q[eldofs]
        state = @view states[:, i]
        state_old = @view states_old[:, i]
        if fielddim>4
            assemble_element_u!(Ke,fe,cellvalues,qe, material,state,state_old)
        else
            assemble_element_phi!(Ke,fe,cellvalues,qe, material,state,state_old)
        end
        assemble!(assembler,eldofs,Ke,fe)
    end
    return K,f
end;

#Subroutine for global assembly - residual only
function assemble_global_r(q::Vector, cellvalues,dh::DofHandler,
                           material::Brittle, states,states_old,store::Bool=true)
    #allocate element stiffness and rhs
    nbase = getnbasefunctions(cellvalues)

    fe = zeros(nbase)

    f=zeros(ndofs(dh))

    #THIS METHOD IS NOT ROBUST
    fielddim=nbase
    for (i, cell) in enumerate(CellIterator(dh))
        reinit!(cellvalues,cell)
        fill!(fe,0)
        eldofs=celldofs(cell)
        qe=q[eldofs]
        state = @view states[:, i]
        state_old = @view states_old[:, i]
        if fielddim>4
            assemble_residual_u!(fe,cellvalues,qe, material,state,state_old,store)
        else
            assemble_residual_phi!(fe,cellvalues,qe, material,state,state_old,store)
        end
        f[eldofs] += fe
    end
    return f
end;

#Compute active strain energy for crack driving force
function CrackDrive(σ::SymmetricTensor{2,3,Float64},ε::SymmetricTensor{2,3,Float64},mat::Brittle)
    flag=mat.flag
    if flag ==Isotropic #Isotropic split
        Psi=0.5*σ ⊡ ε
    elseif flag==VolDev #Volumetric/deviatoric split of Amor
        K = mat.K
        G = mat.G
        Psi = G*dev(ε) ⊡ dev(ε)
        if tr(ε) >0
            Psi += 0.5*K*tr(ε)^2
        end
    elseif flag==Spectral # Spectral split by Miehe
        εₚ= eigvals(ε)
        μ = mat.G
        λ = mat.K - 2μ/3
        Psi = sum(εₚ)>0 ? λ/2*sum(εₚ)^2 : 0.0
        for e in εₚ
            Psi += e>0 ? μ*e^2 : 0.0
        end
    elseif flag== NoTension #Freddis No-tension
        εₚ= sort(eigvals(ε))#note, order is ε₃>ε₂ > ε₁
        μ = mat.G
        λ = mat.K - 2μ/3
        ν = λ/2(λ+μ)

        if εₚ[1] > 0
            Psi = λ/2*sum(εₚ)^2+μ*sum(εₚ.^2)
        elseif ν*εₚ[1]+εₚ[2] > 0
            Psi = λ/2*(εₚ[3]+εₚ[2]+2ν*εₚ[1])^2+μ*((εₚ[3]+ν*εₚ[1])^2+(εₚ[2]+ν*εₚ[1])^2)
        elseif (1-ν)*εₚ[3]+ν*(εₚ[1]+εₚ[2]) > 0
            Psi = λ/2*(1+ν)/(ν*(1-ν^2))*((1-ν)*εₚ[3]+ν*εₚ[1]+ν*εₚ[2])^2
        else
            Psi = 0.0
        end
    end
    return Psi
end

#Track crack growth direction
function CrackTrack(q::Vector, dh::DofHandler, CellValues::CellScalarValues, grid::Grid, a0::Float64, CrackDir::Int64)
    Ac = a0
    v=CrackDir
    for (i, cell) in enumerate(CellIterator(dh))
        reinit!(CellValues,cell)
        eldofs=celldofs(cell)
        ϕe=q[eldofs]
        if maximum(ϕe)>= 0.95
            node_coords=getcoordinates(grid,i)
            for q_point in 1:getnquadpoints(CellValues)
                #Phase field value
                ϕ = function_value(CellValues,q_point,ϕe)
                if ϕ >= 0.95
                    coords = spatial_coordinate(CellValues,q_point,node_coords)
                    Ac = coords[v]>Ac ? coords[v] : Ac
                end
            end
        end
    end

    return Ac
end

#Newton-Raphson solver
function Newton_raphson!(q::Vector, K::SparseMatrixCSC, cellvalues, dh::DofHandler, ch::ConstraintHandler,
                        grid::Grid, Material::Brittle, states, states_old,tag::String)

    K_fac=Cholesky
    iterations = 0
    for nitr = 1:(solver.nitr_inner+1)
        if nitr > solver.nitr_inner
            error("Reached maximum Newton iterations, aborting")
            break
        end
        K, r = assemble_global(q, cellvalues, K, dh, Material,states,states_old);
        apply_zero!(K,r,ch)
        norm_r = maximum(abs.(r[Ferrite.free_dofs(ch)]))
        if (norm_r < solver.NewtonTOL_inner) && (nitr > 1)
            break
        end 
        iterations += 1
        if isposdef(K)
             K_fac = cholesky(K)
        else
            print("matrix not positive definite for $tag \n")
            K_fac = LU(K).L
        end
        Δq = K_fac\r
        q -= Δq
    end
    #    print(tag*" converged in $iterations iterations \n")
    return K_fac, q, iterations
end

#Modified Newton-Raphson solver
function Mod_Newton_Raphson!(q::Vector, K, cellvalues, dh::DofHandler, ch::ConstraintHandler,
                        grid::Grid, Material::Brittle, states, states_old,tag::String)
    iterations = 0
    qold = deepcopy(q)
    fail=false
    for nitr = 1:(solver.nitr_inner+1)
        if nitr > solver.nitr_inner
            fail=true
            q=qold
            break
        end
        r = assemble_global_r(q, cellvalues, dh, Material,states,states_old);
        apply_zero!(r,ch)
        norm_r = maximum(abs.(r[Ferrite.free_dofs(ch)]))
        if (norm_r < solver.NewtonTOL_inner) && (nitr > 1)
            break
        end
        iterations += 1
        Δq = K\r
        q -= Δq
    end
    #    print(tag*" converged in $iterations iterations \n")
    return q, iterations, fail
end

#Determine which solver to use
function DetermineSolver!(solver::SolverState,nitr::Int)
     if nitr==1 && solver.rebuild_start
    
        solver.strategy = NewtonRaphson
        solver.rebuild_counter += 1
        solver.no_restart=0
        solver.rebuild_start = false
        
     elseif solver.rebuild_counter > 3
        solver.strategy = NewtonRaphson
        solver.rebuild_counter += 1
        solver.no_restart=0
        solver.rebuild_start=true
     elseif mod(nitr,solver.rebuild_frequency_base)==0
        solver.strategy = NewtonRaphson
        solver.rebuild_counter += 1
        solver.no_restart=0
     else
        solver.strategy = ModifiedNewtonRaphson
     end
end

function OutputForce(q::Vector, cellvalues, dh::DofHandler,
    grid::Grid, Material::Brittle, states, states_old,set::String)
    F  =assemble_global_r(q, cellvalues, dh, Material,states,states_old, false);
    F_x = 0.0
    F_y = 0.0
    Fnodal = reshape_to_nodes(dh,F,:u)
    if set ∈ keys(grid.nodesets)
        outputset=grid.nodesets[set]
    elseif set ∈ keys(grid.facesets[set])
        print("facesets are currently not supported for force output")
    else
        print("Warning invalid set for force output")
    end
    for  (i,n) in enumerate(grid.nodesets[set])
      F_x += Fnodal[1,n]
      F_y += Fnodal[2,n]
    end
    return F_x, F_y
end

#Main subroutine for processing the finite elements problem.
function Problem(solver::SolverState,output::OutputVariables,cellvalues_ϕ::CellScalarValues, cellvalues_u::CellVectorValues, 
                dh_phi::DofHandler, dh_u::DofHandler,ch_phi::ConstraintHandler, ch_u::ConstraintHandler,
                 grid::Grid, Material::Brittle)

    #MAIN subroutine for controlling the problem solution


    CrackTol = 0.95

    #Create sparsity patterns
    u = zeros(ndofs(dh_u))
    ϕ  =zeros(ndofs(dh_phi))
    K_u=create_sparsity_pattern(dh_u)
    K_ϕ=create_sparsity_pattern(dh_phi)

    #Initialize space for factorized matrices (almost certainly not the best way to do this)
    K_ϕ_fac= Cholesky
    K_u_fac = Cholesky

    # Create material states. One array for each cell, where each element is an array of material-
    # states - one for each integration point
    nqp = getnquadpoints(cellvalues_u)
    states = [MaterialState() for _ in 1:nqp, _ in 1:getncells(grid)]
    states_old = [MaterialState() for _ in 1:nqp, _ in 1:getncells(grid)]


    #This part is necessary for the crack set method
    ∂Ω₃ = getnodeset(grid,"Crack")
    dbc₃ = Dirichlet(:ϕ,∂Ω₃, (x,t) -> 1)
    ch_phi = ConstraintHandler(dh_phi);
    add!(ch_phi,dbc₃)
    close!(ch_phi)
    update!(ch_phi,0.0);
    

    Disp = 0
    #begin computation
    for timestep in 0:n_timesteps
        if solver.loading == Monotonic
            Disp= timestep*solver.u_max/n_timesteps
        elseif solver.loading == Fatigue
            Disp = sin(timestep*3.141592/2.)*solver.u_max
        elseif solver.loading == FatigueCLA
            Disp = timestep== 0 ? 0 : solver.u_max
        end
        update!(ch_u, Disp) # evaluates the D-bndc at time t
        update!(ch_phi)
        apply!(u, ch_u)  # set the prescribed values in the solution vector
        apply!(ϕ, ch_phi)  # set the prescribed values in the solution vector
        solver.rebuild_counter= 0
        solver.no_restart +=1
        #OUTER NR LOOP
        for newton_itr = 1:solver.nitr_outer
            #Too many iterations?
            if newton_itr > solver.nitr_outer
                error("Reached maximum Newton iterations, aborting")
                break
            end
            if newton_itr > 1 #Check convergence
                r_ϕ  =assemble_global_r(ϕ, cellvalues_ϕ, dh_phi, Material,states,states_old);
                norm_r = maximum(abs.(r_ϕ[Ferrite.free_dofs(ch_phi)]))

                if norm_r < solver.NewtonTOL_outer
                    print("\n Time step @time = $timestep, converged in $(newton_itr-1) outer iterations\n")
                    break
                end
            end

            #Decide solution strategy
            DetermineSolver!(solver,newton_itr)
            fail=false
            if solver.strategy == ModifiedNewtonRaphson
                ϕ, nitr_phi, fail = Mod_Newton_Raphson!(ϕ,K_ϕ_fac,cellvalues_ϕ,dh_phi, ch_phi,
                                        grid, Material, states,states_old,"ϕ")
                output.totalIterations_phi += nitr_phi
                
                if !fail
                        u, nitr_u, fail = Mod_Newton_Raphson!(u,K_u_fac,cellvalues_u,dh_u, ch_u,
                                        grid, Material, states,states_old,"u")
                        output.totalIterations_u += nitr_u
                end
                if fail
                    print("Backtracking!")
                    solver.strategy=NewtonRaphson
                end
            end
            if solver.strategy == NewtonRaphson
                print("\n REBUILDING STIFFNESS at time = $timestep\n")
                K_ϕ_fac, ϕ, nitr_phi = Newton_raphson!(ϕ,K_ϕ,cellvalues_ϕ,dh_phi, ch_phi,
                                        grid, Material, states,states_old,"ϕ")
                output.totalIterations_phi += nitr_phi

                K_u_fac, u, nitr_u = Newton_raphson!(u,K_u,cellvalues_u,dh_u, ch_u,
                                        grid, Material, states,states_old,"u")
                output.totalIterations_u += nitr_u
                output.matrixFactorizations_phi += nitr_phi
                output.matrixFactorizations_u += nitr_u
                solver.no_restart=0
                solver.rebuild_counter += 1
            end
            output.totalIterations_outer +=1
        end

        #Evaluate if a rebuild is necessary at beginning of next increment
        if solver.no_restart == solver.rebuild_frequency_min
            solver.rebuild_start=1
        end
        if maximum(ϕ[Ferrite.free_dofs(ch_phi)])≥CrackTol

            #Update Crack Set
            nodesadded=0
            ϕnodal = reshape_to_nodes(dh_phi,ϕ,:ϕ)
            for (i, ϕᵢ) in enumerate(ϕnodal)
                if i ∉ grid.nodesets["Crack"] && ϕᵢ ≥ CrackTol

                    push!(grid.nodesets["Crack"],i)
                    nodesadded +=1
                end
            end
            print("\n Added $nodesadded nodes to crack set at time = $timestep\n")
            ∂Ω₃ = getnodeset(grid,"Crack")
            dbc₃ = Dirichlet(:ϕ,∂Ω₃, (x,t) -> 1.0)
            ch_phi = ConstraintHandler(dh_phi);
            add!(ch_phi,dbc₃)
            close!(ch_phi)
            update!(ch_phi,0.0);
        end  
        #Write output for contour plots
        if mod(timestep,output.plotFrequency)==0 || timestep==n_timesteps
            simtime=timestep/n_timesteps
            saved_file = paraview_collection("TimeSeries"; append = true) do pvd
                vtk_grid("Contour_$(output.plotframe)",dh_phi) do vtk
                    vtk_point_data(vtk,dh_phi,ϕ)
                    vtk_point_data(vtk,dh_u,u)
                    pvd[simtime] = vtk
                end
            end
            output.plotframe += 1
        end
        #Write history output
        if mod(timestep,output.historyFrequency)==0 || timestep==n_timesteps
            #Crack Extension
            CrackExtent= CrackTrack(ϕ, dh_phi,cellvalues_ϕ,grid,output.a0,output.CrackDir)
            F_x, F_y = OutputForce(u, cellvalues_u, dh_u,grid, Material,states,states_old,output.OutputSet)
            historyfile = open("HistoryOutput.txt","a",)
            write(historyfile,"$timestep, $CrackExtent, $(output.matrixFactorizations_u), $(output.matrixFactorizations_phi), "*
                        "$(output.totalIterations_u), $(output.totalIterations_phi), $(output.totalIterations_outer), $Disp, $F_x, $F_y \n")
            close(historyfile)
        end
        states_old .= states
    end
    print("Computation complete! Runtime information:\n")
end;
