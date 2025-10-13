!----------------------------------------------------------------------
! UMAT and UMATHT Subroutines for Phase Field Hydraulic Fracture Modelling
!----------------------------------------------------------------------
!
! This code implements a simple phase field approach for modelling
! hydraulic fracture processes. It includes user subroutines UMAT and UMATHT
! for integration into Abaqus.
!
! License:
! This code is distributed under the BSD license.
!
! Citation:
! If you use this code for academic or industrial purposes, please cite the followings:
!
! Y. Navidtehrani, C. Betegon, E. Martinez-Paneda,
! "A generalised framework for phase field-based modelling of coupled problems:
! Application to thermo-mechanical fracture, hydraulic fracture, hydrogen embrittlement,
! and corrosion," *Engineering Fracture Mechanics*.
! &
! Yousef Navidtehrani, Covadonga Betegón, Javier Vallejos, & Emilio Martínez-Pañeda.
! " A phase field model for hydraulic fracture: Drucker–Prager driving force and a hybrid coupling strategy",
! Computer Methods in Applied Mechanics and Engineering (2025), 444, 118155.
!
! Authors:
! - Yousef Navidtehrani (usofntehrani@gmail.com)
! - Emilio Martinez-Paneda (emilio.martinez-paneda@eng.ox.ac.uk)
!
! Notes:
! - To view mathematical formulas within the code, hover over them using the
!   'Mathover' extension in Visual Studio Code.

c*****************************************************************
!     Transferring variables between subroutines
      module kTransfer

!     For 2D linear elements (CPE4T, CPE3T) and quadratic reduced integration elements (CPE8RT)
      real*8 :: TransE(700000,4) = 0.d0
      real*8 :: TransP(700000,4) = 0.d0
      real*8 :: TransAuxi(700000,4,2) = 0.d0
      real*8 :: TransAlpha(700000,4) = 0.d0

! !     For 2D quadratic full integration elements (CPE8T)
      ! real*8 :: TransE(700000,9) = 0.d0
      ! real*8 :: TransP(700000,9) = 0.d0
      ! real*8 :: TransAuxi(700000,9,2) = 0.d0
      ! real*8 :: TransAlpha(700000,9) = 0.d0

! !     For 3D linear full integration elements (C3D8T)
      ! real*8 :: TransE(700000,8) = 0.d0
      ! real*8 :: TransP(700000,8) = 0.d0
      ! real*8 :: TransAuxi(700000,8,2) = 0.d0
      ! real*8 :: TransAlpha(700000,8) = 0.d0

      real*8 bkk
      
      save
      end module
c*****************************************************************

      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,predef,dpred,
     2cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      use kTransfer
      include 'aba_param.inc'

      character*80 cmname
      character*80 cpname
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),time(2),
     2predef(1),dpred(1),props(nprops),coords(3),drot(3,3),dfgrd0(3,3),
     3dfgrd1(3,3),jstep(4)


!----------------------------------------------------------
! Mechanical and phase field part (MATERIAL-1)
!----------------------------------------------------------  
      if (index(cmname,'MATERIAL-1').eq.1) then

      call getpartinfo(noel, 1, cpname, locnum, jrcd) ! get local element number

!     Initialization
      ddsdde=0.d0
      E=props(1) ! Young's modulus, # Math: E
      xnu=props(2) ! Poisson's ratio, # Math: \nu
      xl=props(3) ! Characteristic length scale, # Math: \ell
      Gc=props(4) ! Fracture energy, # Math: G_c
      kflagS=int(props(5)) ! Solution flag (0: monolithic, 1: staggered)
      AlphaR=props(6) ! Biot’s coefficient of reservoir domain, # Math: \alpha_{\text{r}}
      xc1=props(7) ! First constants for domain indicator fields, # Math: c_1
      xc2=props(8) ! Second constants for domain indicator fields, # Math: c_2
      phi=temp+dtemp ! Phase field variable, # Math: \phi_{t+\delta t}=\phi_t+\delta \phi
      psit=statev(1) ! History variable at the begining of the current increment, # Math: \mathcal{H}
      g=(1.d0-phi)**2+1.d-7 ! Degradation function, # Math: g(\phi)= (1-\phi)^2

!     Build stiffness matrix
      eg=E/(1.d0+xnu)/2.d0 ! Shear modulus, # Math: \mu= \frac{E}{2(1+\nu)}
      elam=(E/(1.d0-2.d0*xnu)-2.d0*eg)/3.d0 ! Lame's first parameter, # Math: \lambda= \frac{E \nu}{(1+\nu)(1-2\nu)}
      bk=E/(1.d0-2.d0*xnu)/3.d0 ! Bulk modulus, # Math: K = \frac{E}{3(1-2\nu)}
      bkk=bk
      do i=1,3
       do j=1,3
        ddsdde(j,i)=elam
       end do
       ddsdde(i,i)=eg*2.d0+elam
      end do
      do i=4,ntens
       ddsdde(i,i)=eg
      end do

!     Update effective stress      
      stran=stran+dstran ! Update strain, # Math: \varepsilon_{t+\delta t} = \varepsilon_t + \delta \varepsilon
      stress=matmul(ddsdde,stran) ! Effective stress of undamaged configuration, # Math: \boldsymbol{\sigma}^{\text{eff}}_0=\boldsymbol{C}_0 \varepsilon

!     Compute the strain energy density
      psi=0.d0
      do i=1,ntens
       psi=psi+stress(i)*stran(i)*0.5d0 ! Strain energy density of undamage configuration, # Math: \psi_0 = \frac{1}{2} \boldsymbol{\sigma_0} : \boldsymbol{\varepsilon}
      end do

!     Apply pressure effect
      !  Indicator fields
      if (phi.le.xc1) then ! Reservoir domain
      chiR=1.d0
      chiF=0.d0
      elseif (phi.ge.xc2) then ! Fracture domain
      chiR=0.d0
      chiF=1.d0
      else ! Transient domain
      chiR=-(phi-xc2)/(xc2-xc1)
      chiF=(phi-xc1)/(xc2-xc1)
      end if
      TransAuxi(locnum,npt,1)=chiR
      TransAuxi(locnum,npt,2)=chiF

      ! Modified Biot's coefficient
      Alpha=chiR*AlphaR+chiF ! # Math: \alpha_b = \chi_{\text{r}} \alpha_{\text{r}} + \chi_{\text{f}}
      TransAlpha(locnum,npt)=Alpha

      stress=stress*g ! Degraded effective stress, # Math: \boldsymbol{\sigma} = g(\phi) \boldsymbol{\sigma}^{\text{eff}}_0
      P=TransP(locnum,npt) ! Pore pressure of current increment from previous iteration, # Math: p_{t+\delta t}^{i-1}
      do i=1,ndi
       stress(i)=stress(i)-Alpha*P ! Update stress, # Math: \boldsymbol{\sigma}= \boldsymbol{\sigma}^{\text{eff}} - \alpha_b p \boldsymbol{I}
      end do

      H=max(psit,psi) ! Phase field history variable, # Math: \mathcal{H} = \max(\psi_0, \mathcal{H})
      ddsdde=ddsdde*g ! Degraded stiffness, # Math: \boldsymbol{C} = g(\phi) \boldsymbol{C}_0
      statev(1)=H
      statev(2)=H
      statev(3)=xl
      statev(4)=Gc

      TransE(locnum,npt)=sum(dstran(1:ndi)) ! Rate of volumetric strain, # Math: \delta {\varepsilon}_{\text{vol}}

      if (kflagS.eq.1) statev(2)=psit ! Staggered scheme

!----------------------------------------------------------
! Fluid part (MATERIAL-2)
!---------------------------------------------------------- 
      elseif (index(cmname,'MATERIAL-2').eq.1) then

      ddsdde=0.d0
      stress=0.d0

      end if

      return
      end

c*****************************************************************

      subroutine umatht(u,dudt,dudg,flux,dfdt,dfdg,
     1statev,temp,dtemp,dtemdx,time,dtime,predef,dpred,
     2cmname,ntgrd,nstatv,props,nprops,coords,pnewdt,
     3noel,npt,layer,kspt,kstep,kinc)

      use kTransfer
      include 'aba_param.inc'

      character*80 cmname
      character*80 cpname
      dimension dudg(ntgrd),flux(ntgrd),dfdt(ntgrd),
     1dfdg(ntgrd,ntgrd),statev(nstatv),dtemdx(ntgrd),
     2time(2),predef(1),dpred(1),props(nprops),coords(3)

!----------------------------------------------------------
! Mechanical and phase field part (MATERIAL-1)
!---------------------------------------------------------- 
      if (index(cmname,'MATERIAL-1').eq.1) then
            
      phi=temp+dtemp ! Phase field variable, # Math: \phi_{t+\delta t}=\phi_t+\delta \phi
      H=statev(2) ! Phase field history variable, # Math: \mathcal{H}
      xl=statev(3) ! Characteristic length scale, # Math: \ell
      Gc=statev(4) ! Fracture energy, # Math: G_c

      U=U+(phi/xl**2-2.d0*(1.d0-phi)*H/(Gc*xl))*DTIME ! Internal energy, # Math: U_{t+\delta t} = U_t + \left(\frac{\phi}{\ell^2} - \frac{2(1-\phi)H}{G_c\ell}\right) \delta t
      DUDT=(1.d0/xl**2+2.d0*H/(Gc*xl))*DTIME ! Derivative of internal energy w.r.t. time, # Math: \frac{\partial U}{\partial \phi} = \left(\frac{1}{\ell^2} + \frac{2H}{G_c\ell}\right) \delta t
      DUDG=0.d0 ! Derivative of internal energy w.r.t. gradient of phase field, # Math: \frac{\partial U}{\partial \nabla \phi} = 0
      DFDT=0.d0 ! Derivative of flux w.r.t. phase field, # Math: \frac{\partial \mathbf{f}}{\partial \phi} = 0
      do i=1,NTGRD
       DFDG(i,i)=-1.d0 ! Derivative of flux w.r.t. gradient of phase field, # Math: \frac{\partial \mathbf{f}}{\partial \nabla \phi} = -\boldsymbol{I}
      end do
      FLUX=matmul(DFDG,DTEMDX) ! Flux, # Math: \mathbf{f} = - \nabla \phi

!----------------------------------------------------------
! Fluid part (MATERIAL-2)
!---------------------------------------------------------- 
      elseif (index(cmname,'MATERIAL-2').eq.1) then

      call getpartinfo(noel, 1, cpname, locnum, jrcd)

!     Initialization
      xRho=props(1) ! Fluid density, # Math: \rho_{\text{fl}}
      xn_pr=props(2) ! Reservoir porosity, # Math: n_{\text{pr}}
      C_fl=props(3) ! Fluid compressibility, # Math: C_{\text{fl}}
      xkR=props(4) ! Reservoir permeability, # Math: {K}_{\text{r}}
      xkF=props(5) ! Fracture permeability, # Math: {K}_{\text{f}}
      xmu=props(6) ! Fluid viscosity, # Math: \mu_{\text{fl}}

!     Updating pressure
      P=temp+dtemp ! Pore pressure at the current increment, # Math: p_{t+\delta t}=p_t+\delta p
      dP=dtemp ! Pore pressure increment, # Math: \delta p

!     Transfer inputs
      dTrE=TransE(locnum,npt) ! Rate of volumetric strain, # Math: \delta {\varepsilon}_{\text{vol}}

!     Indicator fields
      chiR=TransAuxi(locnum,npt,1) ! # Math: {\chi}_{\text{r}}
      chiF=TransAuxi(locnum,npt,2) ! # Math: {\chi}_{\text{f}}

!     Modified Biot's coefficient
      Alpha=TransAlpha(locnum,npt) ! Biot coefficient # Math: \alpha_b

!     Modified porosity
      xn_p=chiR*xn_pr+chiF ! # Math: n_{\text{p}} = \chi_{\text{r}} n_{\text{pr}} + \chi_{\text{f}}

!     Modified permeability tensor
      xKP=(chiR*xkR+chiF*xkF) ! # Math: {K}_{\text{fl}} = \chi_{\text{r}} {K}_{\text{r}} + \chi_{\text{f}} {K}_{\text{f}}

!     Storage coefficient
      StoCoef=(1.d0-Alpha)*(Alpha-xn_p)/bkk+xn_p*C_fl ! # Math: S = \frac{(1 - \alpha_b)(\alpha_b - n_{\text{p}})}{K} + n_{\text{p}} C_{\text{fl}}

!     Fluid equation terms
      U=U+xRho*(StoCoef*dP+Alpha*chiR*dTrE) ! # Math: U_{t+\delta t} =U_{t} +  \left(S \delta p + \alpha_b \chi_{\text{r}} \delta \varepsilon_{\text{vol}}\right)
      DUDT=xRho*StoCoef ! # Math: \frac{\partial U}{\partial p} =S
      DUDG=0.d0 ! # Math: \frac{\partial U}{\partial \nabla p} = 0
      do i=1,NTGRD
      DFDG(i,i)=-xRho*xKP/xmu ! # Math: \frac{\partial \mathbf{f}}{\partial \nabla p} =-\rho_{\mathrm{fl}} \frac{\boldsymbol{K}_{\text{fl}}}{\mu_{\text{fl}}} \boldsymbol{I}
      end do
      DFDT=0.d0 ! # Math: \frac{\partial \mathbf{f}}{\partial p} = 0
      FLUX=matmul(DFDG,DTEMDX) ! # Math: \mathbf{f} =-\rho_{\mathrm{fl}} \frac{\boldsymbol{K}_{\text{fl}}}{\mu_{\text{fl}}} \nabla p

!     Transfering outputs
      TransP(locnum,npt)=P ! Pore pressure, # Math: p
      end if

      return
      end

c*****************************************************************