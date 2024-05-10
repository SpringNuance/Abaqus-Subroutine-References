! This subroutine is trying to replicate the 3 study case results from the paper
! Modelling the coupling between hydrogen diffusion and the mechanical behavior of metals

!***********************************************************************

module common_block
    implicit none
                           
    real*8 :: common_coords(50000, 4, 2)
    real*8 :: grad_sigma_hydrostatic(50000, 4, 2)
    real*8 :: sigma_hydrostatic(50000, 4)

    save 
    ! The save command is very important. 
    ! It allows the values to be stored and shared between subroutines 
    ! without resetting them to zero every time the subroutine is called
end module   

!***********************************************************************

subroutine UEXTERNALDB(lop,lrestart,time,dtime,kstep,kinc)

    use common_block
    include 'aba_param.inc' 
    dimension time(2)
    
    ! LOP=0 indicates that the subroutine is being called at the start of the analysis.
    if (lop == 0) then 
        common_coords = 0.d0
        grad_sigma_hydrostatic = 0.d0
        sigma_hydrostatic = 0.d0
    end if
    
return
end

! CPE8RT: 8-node biquadratic displacement, bilinear temperature, reduced integration

subroutine calculate_grad_sigma_hydrostatic_CPE8RT(noel)

    use common_block

    ! This stores the integration point coordinate in isoparametric spaces
    ! Refers to the progress meeting PPT to know the order of the integration points
    
    dimension deriv(2,4),xjacm(2,2),xjaci(2,2)

    integer, parameter :: isoparametric_coord_X(4) = (/ -1, +1, -1, +1 /)
    integer, parameter :: isoparametric_coord_Y(4) = (/ -1, -1, +1, +1 /)
    
    do int_id = 1,4 
        ! print *, 'int_id: ', int_id
        s = isoparametric_coord_X(int_id)
        t = isoparametric_coord_Y(int_id)
        
        deriv(1,1) = -(1.d0/4.0) * (1-t)
        deriv(1,2) =  (1.d0/4.0) * (1-t)
        deriv(1,3) = -(1.d0/4.0) * (1+t)
        deriv(1,4) =  (1.d0/4.0) * (1+t)
        deriv(2,1) = -(1.d0/4.0) * (1-s)
        deriv(2,2) = -(1.d0/4.0) * (1+s)
        deriv(2,3) =  (1.d0/4.0) * (1-s)
        deriv(2,4) =  (1.d0/4.0) * (1+s)

        xjacm(1,1) = deriv(1,1) * common_coords(noel,1,1) + deriv(1,2) * common_coords(noel,2,1) &
                   + deriv(1,3) * common_coords(noel,3,1) + deriv(1,4) * common_coords(noel,4,1)
    
        xjacm(1,2) = deriv(1,1) * common_coords(noel,1,2) + deriv(1,2) * common_coords(noel,2,2) &
                   + deriv(1,3) * common_coords(noel,3,2) + deriv(1,4) * common_coords(noel,4,2)
    
        xjacm(2,1) = deriv(2,1) * common_coords(noel,1,1) + deriv(2,2) * common_coords(noel,2,1) &
                   + deriv(2,3) * common_coords(noel,3,1) + deriv(2,4) * common_coords(noel,4,1)
    
        xjacm(2,2) = deriv(2,1) * common_coords(noel,1,2) + deriv(2,2) * common_coords(noel,2,2) &
                   + deriv(2,3) * common_coords(noel,3,2) + deriv(2,4) * common_coords(noel,4,2)

        djacb = xjacm(1,1) * xjacm(2,2) - xjacm(1,2) * xjacm(2,1) 
    
        xjaci(1,1) =  xjacm(2,2)/djacb 
        xjaci(1,2) = -xjacm(1,2)/djacb  
        xjaci(2,1) = -xjacm(2,1)/djacb   
        xjaci(2,2) =  xjacm(1,1)/djacb

        a1 = xjaci(1,1) * deriv(1,1) + xjaci(1,2) * deriv(2,1) 
        a2 = xjaci(1,1) * deriv(1,2) + xjaci(1,2) * deriv(2,2) 
        a3 = xjaci(1,1) * deriv(1,3) + xjaci(1,2) * deriv(2,3) 
        a4 = xjaci(1,1) * deriv(1,4) + xjaci(1,2) * deriv(2,4) 
        b1 = xjaci(2,1) * deriv(1,1) + xjaci(2,2) * deriv(2,1) 
        b2 = xjaci(2,1) * deriv(1,2) + xjaci(2,2) * deriv(2,2) 
        b3 = xjaci(2,1) * deriv(1,3) + xjaci(2,2) * deriv(2,3) 
        b4 = xjaci(2,1) * deriv(1,4) + xjaci(2,2) * deriv(2,4)  
    
        grad_sigma_hydrostatic(noel,int_id, 1) = a1 * sigma_hydrostatic(noel,1) + &
                                                 a2 * sigma_hydrostatic(noel,2) + &
                                                 a3 * sigma_hydrostatic(noel,3) + &
                                                 a4 * sigma_hydrostatic(noel,4)
        grad_sigma_hydrostatic(noel,int_id, 2) = b1 * sigma_hydrostatic(noel,1) + &
                                                 b2 * sigma_hydrostatic(noel,2) + &
                                                 b3 * sigma_hydrostatic(noel,3) + &
                                                 b4 * sigma_hydrostatic(noel,4)
    end do 

return
end

!***********************************************************************

! This is isotropic von Mises plasticity model
subroutine UMAT(stress,statev,ddsdde,sse,spd,scd, &
    rpl,ddsddt,drplde,drpldt, &
    stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname, &
    ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt, &
    celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc) 
!
    use common_block
    include 'aba_param.inc'
!
    character*80 cmname
    dimension stress(ntens),statev(nstatv), &
       ddsdde(ntens,ntens), &
       ddsddt(ntens),drplde(ntens), &
       stran(ntens),dstran(ntens),time(2),predef(1),dpred(1), &
       props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
!
    real E, nu, sigma_0, n_hardening, lambda, mu, PEEQ, von_Mises_stress, hydrostatic_stress 
    real sigma_Y, sigma_H0, dPEEQ, E_tangent, effective_mu, effective_lambda, effective_hard

    dimension eps_elastic(ntens),eps_plastic(ntens),flow_stress(ntens), &
         old_stress(ntens), old_eps_platic(ntens)


    real, parameter :: toler = 1.d-6
    real, parameter :: newton = 50

    ! Initialize material properties
    
    E = props(1)           ! Young's modulus
    nu = props(2)          ! Poisson's ratio
    sigma_0 = props(3)     ! Initial yield strength in the absence of hydrogen (MPa)
    n_hardening = props(4) ! Strain hardening exponent not affected by Hydrogen
                           ! in paper it is 5, but it is already inversed here, so it is 0.2
    dPEEQ = 0.d0           ! Equivalent plastic strain increment


    old_stress = stress
    old_eps_platic = eps_plastic

! UMAT and UMATHT are integration point level, 
! which loops over all the integration points over all elements

! noel: The current element number
! npt: The current integration point number

! TIME(1)
! Value of step time at the beginning of the current increment or frequency.

! TIME(2)
! Value of total time at the beginning of the current increment.

!  Compute the gradient of the hydrostatic stress    

    if (npt == 1 .and. time(1) > 0) then
        call calculate_grad_sigma_hydrostatic_CPE8RT(noel)
    end if     
      
! Lame's parameters
    mu = E/(2.0d0 * (1.0 + nu))  ! Shear modulus
    !lambda = E*nu/((1.0 + nu) * (1.0 - 2.0 * nu)) ! Lame's constant
    lambda = (E/(1.d0-2.d0*nu)-2.d0*mu)/3.d0 ! Lame's constant
! Stiffness matrix
    
    ! initialize as 0
    ddsdde = 0.0
    
    do i = 1, ndi
        do j = 1, ndi
            ddsdde(j, i) = lambda
        end do 
        ddsdde(i,i) = lambda + 2.0 * mu
    end do 

! Shear contribution
    do i = ndi + 1, ntens
        ddsdde(i,i) = mu
    end do 

    call rotsig(statev(1), drot, eps_elastic, 2, ndi, nshr)
    call rotsig(statev(ntens+1), drot, eps_plastic, 2, ndi, nshr)
    
    PEEQ = statev(9)
    
! Stress increment evaluation
    stress = stress + matmul(ddsdde,dstran) 

! calculate elastic strain
    eps_elastic = eps_elastic + dstran

! Calculate equivalent von Mises stress
    
    von_Mises_stress = (stress(1) - stress(2))**2 + &
                       (stress(2) - stress(3))**2 + &
                       (stress(3) - stress(1))**2
    do i = ndi + 1, ntens
        von_Mises_stress = von_Mises_stress + 6.d0 * stress(i)**2
    end do
    von_Mises_stress = sqrt(von_Mises_stress/2.d0)

    ! Extracting PSI_Cbar_L from the statev array
    ! PSI_Cbar_L = statev(13)

    ! First extreme case: PSI_Cbar_L = 1.0
    ! PSI_Cbar_L = 1.0
    
    ! Second extreme case: PSI_Cbar_L = 0.2 (xi)
    ! PSI_Cbar_L = 0.2
    
    if (time(1) > 0) then
        PSI_Cbar_L = statev(13)
    else
        PSI_Cbar_L = 0.52
    end if

    sigma_H0 = PSI_Cbar_L * sigma_0 ! constant within this UMAT
    
    ! From now, sigma_H0 is the yield strength in the presence of hydrogen

    ! Equation 33: We assume that the current
    ! yield strength, sigma_Y, is a function of the equivalent plastic strain PEEQ
    ! and Cbar_L according to the relationship
  
    sigma_Y = sigma_H0 * ((1 + E * PEEQ/sigma_H0) ** n_hardening)

    ! Determine if active yielding
    if (von_Mises_stress > (1.d0 + toler) * sigma_Y) then

        ! Calculate the flow_stress direction
        hydrostatic_stress = (stress(1) + stress(2) + stress(3))/3.d0
        flow_stress(1:ndi) = (stress(1:ndi) - hydrostatic_stress)/von_Mises_stress
        flow_stress(ndi+1:ntens) = stress(ndi+1:ntens)/von_Mises_stress
        
        ! Newton-Raphson iterative solution
        ! the Newton-Raphson method aims to solve the following nonlinear equation:
        ! von_Mises_stress - (3.d0 * mu * dPEEQ) - sigma_Y = 0
        ! Newton-Raphson is an iterative root-finding method used to solve equations of the form 
        ! f(x)=0. It uses the derivative of the function f(x) to iteratively converge to a root. 
        ! The method requires an initial guess and proceeds as follows:

        ! x_{n+1} = x_n - f(x_n)/f'(x_n)
        ! where x_n is the current guess and x_{n+1} is the next guess.
        ! f(x_n) is the function evaluated at the current estimate.
        ! f'(x_n) is the derivative of the function evaluated at the current estimate.

        ! initial guess for dPEEQ
        dPEEQ = 0.d0
    
        ! Tangent Modulus Derivation
        ! It represents the slope of the stress-strain curve at any point, particularly in the plastic region of the curve.
        ! It is the derivative of the yield stress with respect to the plastic strain (derivative of sigma_Y w.r.t PEEQ)
        E_tangent = E * n_hardening * (1.d0 + E * PEEQ/sigma_H0) ** (n_hardening - 1)
        
        do k_newton = 1, newton
            rhs = von_Mises_stress - (3.d0 * mu * dPEEQ) - sigma_Y
            dPEEQ = dPEEQ + rhs / ((3.d0 * mu) + E_tangent)

            sigma_Y = sigma_H0 * (1.d0 + E * (PEEQ + dPEEQ)/sigma_H0) ** n_hardening
            E_tangent = E * n_hardening * (1.d0 + E * (PEEQ + dPEEQ)/sigma_H0) ** (n_hardening-1)
            
            if (abs(rhs) < toler * sigma_H0) exit
        end do

        if (k_newton == newton) write(7,*)'WARNING: plasticity loop failed'

        ! Update stresses and strains
        stress(1:ndi) = flow_stress(1:ndi) * sigma_Y + hydrostatic_stress

        ! Update the elastic and plastic strains
        eps_plastic(1:ndi) = eps_plastic(1:ndi) + 3.d0/2.d0 * flow_stress(1:ndi) * dPEEQ
        eps_elastic(1:ndi) = eps_elastic(1:ndi) - 3.d0/2.d0 * flow_stress(1:ndi) * dPEEQ
        
        stress(ndi+1:ntens) = flow_stress(ndi+1:ntens) * sigma_Y

        eps_plastic(ndi+1:ntens) = eps_plastic(ndi+1:ntens) + 3.d0 * flow_stress(ndi+1:ntens) * dPEEQ
        eps_elastic(ndi+1:ntens) = eps_elastic(ndi+1:ntens) - 3.d0 * flow_stress(ndi+1:ntens) * dPEEQ

        ! Finally, we update the equivalent plastic strain
        PEEQ = PEEQ + dPEEQ

!       Calculate the plastic strain energy density
        do i = 1, ntens
            spd = spd+(stress(i) + old_stress(i)) * (eps_plastic(i) - old_eps_platic(i))/2.d0
        end do

!       Formulate the jacobian (material tangent)   

        ! effective shear modulus
        effective_mu = mu * sigma_Y / von_Mises_stress 

        ! effective Lame's constant
        effective_lambda = (E/(1.d0 - 2.d0 * nu) - 2.d0 * effective_mu)/3.d0 

        ! effective hardening modulus
        effective_hard = 3.d0 * mu * E_tangent/(3.d0 * mu + E_tangent) - 3.d0 * effective_mu 

        do i = 1, ndi
            do j = 1, ndi
                ddsdde(j,i) = effective_lambda
            end do
            ddsdde(i,i) = 2.d0 * effective_mu + effective_lambda
        end do

        do i = ndi + 1, ntens
            ddsdde(i,i) = effective_mu
        end do

        do i = 1, ntens
            do j = 1, ntens
                ddsdde(j,i) = ddsdde(j,i) + effective_hard * flow_stress(j) * flow_stress(i)
            end do
        end do
    endif


!   Store strains in state variable array (ntens = 4 in our case)
    statev(1:4) = eps_elastic           
    statev(5:8) = eps_plastic          
    statev(9) = PEEQ                    
    statev(10) = dPEEQ                               
    statev(11) = (stress(1) + stress(2) + stress(3))/3.d0 ! Hydrostatic stress
    
! Storing the hydrostatic stress and the common_coords of the current integration point   

! npt is from 1 to 9
    sigma_hydrostatic(noel, npt) = (stress(1) + stress(2) + stress(3)) / 3.0      
    common_coords(noel, npt, 1) = coords(1)
    common_coords(noel, npt, 2) = coords(2)

return
end

!***********************************************************************

subroutine UMATHT(u,dudt,dudg,flux,dfdt,dfdg, &
    statev,temp,dtemp,dtemdx,time,dtime,predef,dpred, &
    cmname,ntgrd,nstatv,props,nprops,coords,pnewdt, &
    noel,npt,layer,kspt,kstep,kinc)
! 
    ! include 'common_block.inc'
    use common_block
    ! use sleep_module
    include 'aba_param.inc'
!
    character(len=80) :: cmname
    dimension dudg(ntgrd),flux(ntgrd),dfdt(ntgrd), &
      dfdg(ntgrd,ntgrd),statev(nstatv),dtemdx(ntgrd), &
      time(2),predef(1),dpred(1),props(nprops),coords(3)
    
! user coding to define u,dudt,dudg,flux,dfdt,dfdg,
! and possibly update statev, pnewdt
    
    ! Define all double precision for all variables used 
    double precision :: K ! Arhenius reaction rate constant

    double precision :: NA = 6.023D23 ! Avogadro's number (1/mol)
    double precision :: N_L = 9.24D28 ! Number of solvent atoms (Ni) per unit volume (1/m^3)
    double precision :: beta = 6.0D0 ! Number of interstitial sites per solvent (Ni) atom
    
    double precision :: a_lattice = 2.86D-10 ! Lattice parameter (m)
    double precision :: WB = -18.0D3 ! Binding energy of hydrogen to dislocations (J/mol)

    double precision :: xi = 0.2D0 ! xi fitting parameter in Eq 37 

    double precision :: Cbar_max = 35.0D0 ! (mol/m^3) in Eq 37
    double precision :: Cbar_min = 15.0D0 ! (mol/m^3) in Eq 37

    double precision :: gamma = 2.0D16 ! gamma fitting parameter in Eq 35 (1/m^2)
    double precision :: rho_d0 = 10.0D10 ! Dislocation density for the annealed material (1/m^2)

    R        = props(1) ! Universal gas constant J/(mol K)
    T        = props(2) ! Temperature (K)
    VH       = props(3) ! part molar volume of hydrogen (m^3/mol), 
    DL       = props(4) ! Diffusion coefficient for hydrogen in Nickel (m^2/s)   
     
    ! Since dudt is dCbar_total/dCbar_L, current temperature is the current Cbar_L
    ! and the current dtemp is the current dCbar_L 
    
    ! Obtain the current Cbar_L
    ! prev_Cbar_L = temp
    ! dCbar_L = dtemp
    Cbar_L = temp + dtemp

    ! Equation 15: Arhenius reaction rate constant 
    K = exp( -WB / (R * T)) ! constant around 1361.5 (dimless)

    ! Finding theta_L based on equations (1) and (5)
    C_L = Cbar_L * NA
    theta_L = C_L / (beta * N_L) 

    ! Finding theta_trap based on Oriani equilibrium theory, which results in a Fermi-Dirac relation
    ! theta_trap / (1 - theta_trap) = K * theta_L / (1 - theta_L)
    ! However if theta_L << 1  then 
    ! theta_trap / (1 - theta_trap) = K * theta_L (Equation 14)
    temp_result = K * theta_L
    theta_trap = temp_result / (1 + temp_result)
    
    ! Extract PEEQ and calculate rho_d and dNbar_trap_dPEEQ
    PEEQ = statev(9)
    
    if (PEEQ < 0.5) then
        ! Equation 35: Hydrogen trapping in dislocations
        rho_d = rho_d0 + gamma * PEEQ
        ! Equation 4: We assume that N_trap is proportional to the dislocation density rho_d
        ! and that the dislocation density is a function of the accumulated plastic strain PEEQ
        N_trap = (sqrt(2.0)/lattice_a) * rho_d 
        Nbar_trap = N_trap / NA
        ! Equation 17: relationship between Nbar and PEEQ
        ! = 1.414 / 2.86 * 10^(-10) * 2 * 10^16 / 6.023 * 10^23 = 164.197
        dNbar_trap_dPEEQ = gamma * (sqrt(2.0) / a_lattice) / NA 
    else
        rho_d = 10D16 ! Equation 36
        N_trap = (sqrt(2.0)/lattice_a) * rho_d
        Nbar_trap = N_trap / NA
        dNbar_trap_dPEEQ = 0
    end if

    ! Equation 17: obtain Nbar_L
    Nbar_L = N_L / NA
    
    ! Equation 4: obtain C_trap
    C_trap = 1 * theta_trap * N_trap

    ! Equation 6: obtain Cbar_trap
    Cbar_trap = C_trap / NA

    ! Equation (23) The 1st equation to update dudt (partial_Cbar_total_partial_Cbar_L)
    dudt = 1 + (beta * Nbar_trap * K * Nbar_L) / ((beta * Nbar_L + K * Cbar_L) ** 2)
    
    ! Equation (23) The 2nd equation to update dudtrap (partial_Cbar_total_partial_Nbar_trap)
    dudtrap = (K * Cbar_L) / (K * Cbar_L + beta * Nbar_L)

    ! Equation (22) The total hydrogen diffusion equation to update u 
    ! Cbar_total_t+1 = Cbar_total_t + partial_Cbar_total_partial_Cbar_L * dCbar_L 
    !                               + partial_Cbar_total_partial_Nbar_trap * dNbar_trap_dPEEQ * dPEEQ
    dPEEQ = statev(10) 
    
    u = u + dudt * dtemp + dudtrap * dNbar_trap_dPEEQ * dPEEQ
    Cbar_total = u
    
    ! Since the problem is 2-dimensional, ntgrd = 2
    ! ntgrd: Number of spatial gradients of temperature

    do i = 1, ntgrd

        ! Equation (10) to update the flux
        grad_Cbar_L_i = dtemdx(i)
        flux(i) = ((DL * Vh * Cbar_L) / (R * T)) * grad_sigma_hydrostatic(noel, npt, i) &
                                                 - DL * grad_Cbar_L_i 
        ! Assumed to be 0
        dudg(i) = 0.0

        ! Equation (23) The 3rd equation to update dfdt
        ! part_Jm_part_Cbar_L = (DL * VH) / (R * T) * grad_sigma_hydrostatic(noel, npt, i) ! from common block
        dfdt(i) = (DL * VH) / (R * T) * grad_sigma_hydrostatic(noel, npt, i)

        ! Equation (23) The 4th equation to update dfdg
        ! part_Jm_part_grad_Cbar_L = - DL
        dfdg(i,i) = -DL

    end do
    
    ! Updating the state variables
    
    ! Equation 37: monotonically decreasing function of Cbar_L, which will then be passed to UMAT
    PSI_Cbar_L = 1 - (1 - xi) * ((Cbar_L - Cbar_min) / (Cbar_max - Cbar_min))
    
    statev(12) = rho_d
    statev(13) = PSI_Cbar_L
    statev(14) = Cbar_L
    statev(15) = Cbar_trap
    statev(16) = Cbar_total
    statev(17) = theta_L
    statev(18) = theta_trap

    return
    end