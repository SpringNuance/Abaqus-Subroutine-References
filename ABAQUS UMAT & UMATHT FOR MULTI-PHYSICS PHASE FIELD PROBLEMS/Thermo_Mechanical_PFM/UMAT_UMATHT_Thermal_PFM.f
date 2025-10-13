!----------------------------------------------------------------------
! UMAT and UMATHT Subroutines for Phase Field Thermal Fracture Modelling
!----------------------------------------------------------------------
!
! This code implements a simple phase field approach for modelling
! thermo-mechanical fracture processes. It includes user subroutines UMAT and UMATHT
! for integration into Abaqus.
!
! License:
! This code is distributed under the BSD license.
!
! Citation:
! If you use this code for academic or industrial purposes, please cite the following:
!
! Y. Navidtehrani, C. Betegon, E. Martinez-Paneda,
! "A generalised framework for phase field-based modelling of coupled problems:
! Application to thermo-mechanical fracture, hydraulic fracture, hydrogen embrittlement,
! and corrosion," *Engineering Fracture Mechanics*.
!
! Authors:
! - Yousef Navidtehrani (usofntehrani@gmail.com)
! - Emilio Martinez-Paneda (emilio.martinez-paneda@eng.ox.ac.uk)
!
! Notes:
! - To view mathematical formulas within the code, hover over them using the
!   'Mathover' extension in Visual Studio Code.

c**********************************************************************
!     Transfering variables between subroutines
      module kTransfer

!     For 2D linear elements (CPE4T, CPE3T) and quadratic reduced integration elements (CPE8RT)
       real*8 :: Transg(700000,4) = 0.d0
       real*8 :: TransT_0(700000,4) = 0.d0
       real*8 :: TransT(700000,4) = 0.d0

! !     For 2D quadratic full integration elements (CPE8T)
!        real*8 :: Transg(700000,9) = 0.d0
!        real*8 :: TransT_0(700000,9) = 0.d0
!        real*8 :: TransT(700000,9) = 0.d0

! !     For 3D linear full integration elements (C3D8T)
!        real*8 :: Transg(700000,8) = 0.d0
!        real*8 :: TransT_0(700000,8) = 0.d0
!        real*8 :: TransT(700000,8) = 0.d0
       
      save
      end module
c**********************************************************************

      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,predef,dpred,
     2cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      use ktransfer
      include 'aba_param.inc'

      character*80 cmname
      character*80 cpname
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),time(2),
     2predef(1),dpred(1),props(nprops),coords(3),drot(3,3),dfgrd0(3,3),
     3dfgrd1(3,3),jstep(4)

!----------------------------------------------------------
! Mechanical and phase field part (MAT-MECH)
!----------------------------------------------------------    
      if (index(cmname,'MAT-MECH').eq.1) then

      call getpartinfo(noel, 1, cpname, locnum, jrcd) ! get local element number

!     Initialization
      ddsdde=0.d0
      E=props(1) ! Young's modulus, # Math: E
      xnu=props(2) ! Poisson's ratio, # Math: \nu
      alpha_T=props(3) ! Thermal expansion coefficient, # Math: \alpha_T
      kflagS=int(props(4)) ! Solution flag (0: monolithic, 1: staggered)
      phi=temp+dtemp ! Phase field variable, # Math: \phi
      psit=statev(1) ! History variable at the begining of the current increment, # Math: \mathcal{H}
      g=(1.d0-phi)**2 ! Degradation function, # Math: g(\phi)= (1-\phi)^2

!     Import field variables from other subroutines
      T_0=TransT_0(locnum,npt)! Initial temperature, # Math: T_{0} 
      T=TransT(locnum,npt)! Temperature from previous iteration, # Math: T

!     Build stiffness matrix
      eg=E/(1.d0+xnu)/2.d0 ! Shear modulus, # Math: \mu= \frac{E}{2(1+\nu)}
      elam=(E/(1.d0-2.d0*xnu)-2.d0*eg)/3.d0 ! Lame's first parameter, # Math: \lambda= \frac{E \nu}{(1+\nu)(1-2\nu)}
      do i=1,3
       do j=1,3
        ddsdde(j,i)=elam
       end do
       ddsdde(i,i)=eg*2.d0+elam
      end do
      do i=4,ntens
       ddsdde(i,i)=eg
      end do

!     Update stresses      
      stran=stran+dstran ! Update strain, # Math: \varepsilon_{t+\delta t} = \varepsilon_t + \delta \varepsilon
      do i=1,ndi
        stran(i)=stran(i)-alpha_T*(T-T_0) ! Elastic strain, # Math: \varepsilon_{e} =\varepsilon-\varepsilon_T=\varepsilon -\alpha_T (T-T_0) \boldsymbol{I}
      end do
      stress=matmul(ddsdde,stran) ! Undamaged stress, # Math: \sigma_0 = C \varepsilon_e

!     Compute the strain energy density
      psi=0.d0
      do i=1,ntens
       psi=psi+stress(i)*stran(i)*0.5d0 ! Strain energy density of undamage configuration, # Math: \psi_0 = \frac{1}{2} \sigma_0 : \varepsilon
      end do
      H=max(psit,psi) ! Phase field history variable, # Math: \mathcal{H} = \max(\psi_0, \mathcal{H})
      stress=stress*g ! Degraded stress, # Math: \sigma = g(\phi) \sigma_0
      ddsdde=ddsdde*g ! Degraded stiffness, # Math: C = g(\phi) C_0
      statev(1)=H 
      statev(2)=H

      if (kflagS.eq.1) statev(2)=psit ! Staggered scheme

!     Export field variables for other subroutines
      Transg(locnum,npt)=g ! Degradation function

!----------------------------------------------------------
! Thermal part (MAT-HEAT)
!---------------------------------------------------------- 
      else if (index(cmname,'MAT-HEAT').eq.1) then

      stress=0.d0
      ddsdde=0.d0
     
      end if

      return
      end

c**********************************************************************

      subroutine umatht(u,dudt,dudg,flux,dfdt,dfdg,
     1statev,temp,dtemp,dtemdx,time,dtime,predef,dpred,
     2cmname,ntgrd,nstatv,props,nprops,coords,pnewdt,
     3noel,npt,layer,kspt,kstep,kinc)

      use ktransfer 
      include 'aba_param.inc'

      character*80 cmname
      character*80 cpname
      dimension dudg(ntgrd),flux(ntgrd),dfdt(ntgrd),
     1dfdg(ntgrd,ntgrd),statev(nstatv),dtemdx(ntgrd),
     2time(2),predef(1),dpred(1),props(nprops),coords(3)

!----------------------------------------------------------
! Mechanical and phase field part (MAT-MECH)
!---------------------------------------------------------- 
      if (index(cmname,'MAT-MECH').eq.1) then

      call getpartinfo(noel, 1, cpname, locnum, jrcd) ! get local element number

      xl=props(1) ! Characteristic length scale, # Math: \ell
      Gc=props(2) ! Fracture energy, # Math: G_c
      phi=temp+dtemp ! Phase field variable, # Math: \phi
      H=statev(2) ! Phase field history variable, # Math: \mathcal{H}

      U=U+(phi/xl**2-2.d0*(1.d0-phi)*H/(Gc*xl))*dtime ! Internal energy, # Math: U_{t+\delta t} = U_t + \left(\frac{\phi}{\ell^2} - \frac{2(1-\phi)H}{G_c\ell}\right) \delta t
      DUDt=(1.d0/xl**2+2.d0*H/(Gc*xl))*dtime ! Derivative of internal energy w.r.t. time, # Math: \frac{\partial U}{\partial t} = \left(\frac{1}{\ell^2} + \frac{2H}{G_c\ell}\right) \delta t
      DUDg=0.d0 ! Derivative of internal energy w.r.t. gradient of phase field, # Math: \frac{\partial U}{\partial \nabla \phi} = 0
      DFDT=0.d0 ! Derivative of flux w.r.t. phase field, # Math: \frac{\partial \mathbf{f}}{\partial \phi} = 0
      do i=1,ntgrd
       DFDG(i,i)=-1.d0 ! Derivative of flux w.r.t. gradient of phase field, # Math: \frac{\partial \mathbf{f}}{\partial \nabla \phi} = -\boldsymbol{I}
      end do
      FLUX=matmul(DFDG,DTEMDX) ! Flux, # Math: \mathbf{f} = - \nabla \phi

!----------------------------------------------------------
! Thermal part (MAT-HEAT)
!---------------------------------------------------------- 
      else if (index(cmname,'MAT-HEAT').eq.1) then

      call getpartinfo(noel, 1, cpname, locnum, jrcd) ! get local element number

      xk_T=props(1) ! Thermal conductivity, # Math: k_T
      c_T=props(2) ! Specific heat, # Math: c_T
      kflagT=int(props(3)) ! Method of applying degredation function (1: non-degraded thermal conductivity, 2: degraded thermal conductivity)
      if ((kflagT.ne.1).and.(kflagT.ne.2)) then
        write(6,*) 'Error: props(3) (kflagT) must be either 1 or 2.'
        call xit
      end if ! Check if thermal conductivity is degraded

!     Import field variables from other subroutines
      g=Transg(locnum,npt) ! Degradation function

      if (kflagT.eq.2) xk_T=g*xk_T ! Degraded thermal conductivity, # Math: k_T = g(\phi) k_{0}

      U=U+(c_T*dtemp) ! Internal energy, # Math: U_{t+\delta t} = U_t + c_T \delta T
      DUDT=c_T ! Derivative of internal energy w.r.t. time, # Math: \frac{\partial U}{\partial T} = c_T
      DUDG=0.d0 ! Derivative of internal energy w.r.t. gradient of phase field, # Math: \frac{\partial U}{\partial \nabla T} = 0
      DFDT=0.d0 ! Derivative of flux w.r.t. phase field, # Math: \frac{\partial \mathbf{f}}{\partial T} = 0
      do i=1,ntgrd
      DFDG(i,i)=-xk_T ! Derivative of flux w.r.t. gradient of phase field, # Math: \frac{\partial \mathbf{f}}{\partial \nabla T} = xk_T \boldsymbol{I}
      end do
      FLUX=matmul(DFDG,DTEMDX) ! Flux, # Math: \mathbf{f} = xk_T \nabla T

!     Export field variables for other subroutines
      if ((time(2)).eq.0.d0) TransT_0(locnum,npt)=temp ! Initial temperature, # Math: T_{0}
      TransT(locnum,npt)=temp+dtemp ! Current temperature, # Math: T_{t+\delta t} = T_t + \delta T 

      end if

      return
      end

c**********************************************************************