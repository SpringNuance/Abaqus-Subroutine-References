!----------------------------------------------------------------------
! UMAT and UMATHT Subroutines for Phase Field Hydrogen Embrittlement Modelling
!----------------------------------------------------------------------
!
! This code implements a simple phase field approach for modelling
! hydrogen embrittlement processes. It includes user subroutines UMAT and UMATHT
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

      module ktransfer

      real*8 :: SH(500000,9) = 0.d0 
      real*8 :: GSH(500000,9,3)=0.d0 
      real*8 :: transcH(500000,9)= 0.d0
      real*8 :: CoordNode(500000,2)=0.d0
      real*8 :: TransBx(500000,9,8)=0.d0
      real*8 :: TransBy(500000,9,8)=0.d0

      integer :: ElemConnec(500000,8)=0
      integer numElem, kflagS

      save
      end module

c*****************************************************************

      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)

      use ktransfer
      INCLUDE 'ABA_PARAM.INC'

      DIMENSION TIME(2)
      character*256 JOBNAME

      dimension dNg_1(8,8), Var(1,8), GVar(2)

!     Computing shape function and it's derivative at the start of the analysis
      if (lop.eq.0) then
            CALL GETJOBNAME( JOBNAME, LENJOBNAME )
            call KNodeElemInfo(JOBNAME)
            call kInvNg(dNg_1)
            do noel=1,numElem
             do npt=1,9
              call kjacobian(noel,npt,dNg_1)
             end do
            end do
          
!     At the start of the current analysis increment
      elseif (lop.eq.1) then

      if (kflagS.eq.1) then
       do locnum=1,numElem

!     Computing gradient of damaged hydrostatic stress
        Var(1,1:8)=SH(locnum,1:8)
        do npt=1,9
         call kGardVar(locnum,npt,Var,GVar)
         do i=1,2
          GSH(locnum,npt,i)=GVar(i)
         end do
        end do
      end do
      end if

      end if

      RETURN
      END

c*****************************************************************

      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,predef,dpred,
     2cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      use ktransfer
      include 'aba_param.inc'

      character*80 cmname
      CHARACTER*80 CPNAME
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),time(2),
     2predef(1),dpred(1),props(nprops),coords(3),drot(3,3),dfgrd0(3,3),
     3dfgrd1(3,3),jstep(4)

      dimension Var(1,8), GVar(2)


!!!!!! Part 1 including mechanical and phase field DOFs
      if (cmname.eq.'MATERIAL-1') then

!     Finding local element number
       CALL GETPARTINFO(noel, 1, CPNAME, locnum, JRCD)

!     Initialization
       ddsdde=0.d0
       E=props(1) ! Young's modulus, # Math: E
       xnu=props(2) ! Poisson's ratio, # Math: \nu
       kflagS=int(props(3)) ! Solution flag (0: monolithic, 1: staggered)
       phi=temp+dtemp ! Phase field variable, # Math: \phi_{t+\delta t}=\phi_t+\delta \phi
       psit=statev(1) ! History variable at the beginning of the current increment, # Math: \mathcal{H}
       g=(1.d0-phi)**2+1.d-07 ! Degradation function, # Math: g(\phi)= (1-\phi)^2

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
       stran=stran+dstran ! Update strain, # Math: \boldsymbol{\varepsilon_{t+\delta t} = \boldsymbol{\varepsilon}_t + \delta \boldsymbol{\varepsilon}}
       stress=matmul(ddsdde,stran) ! Stress of undamaged configuration, # Math: \boldsymbol{\sigma}_0=\boldsymbol{C}_0 \boldsymbol{\varepsilon}

!     Compute the strain energy density
       psi=0.d0
       do i=1,ntens
        psi=psi+stress(i)*stran(i)*0.5d0 ! Strain energy density of undamage configuration, # Math: \psi_0 = \frac{1}{2} \boldsymbol{\sigma_0} : \boldsymbol{\varepsilon}
       end do
       H=max(psit,psi) ! Phase field history variable, # Math: \mathcal{H} = \max(\psi_0, \mathcal{H})

!     damaged hydrostatic stress
       SH(locnum,npt)=(1.d0-phi)**2*(stress(1)+stress(2)+stress(3))/3.d0 ! Degraded hydrostatic stress, # Math: \sigma_h=g(\phi) \frac{\sigma_{11}+\sigma_{22}+\sigma_{33}}{3}
!     Applying degradation into stress and stiffness
       stress=stress*g ! Degraded stress, # Math: \boldsymbol{\sigma} = g(\phi) \boldsymbol{\sigma}_0
       ddsdde=ddsdde*g ! Degraded stiffness, # Math: \boldsymbol{C} = g(\phi) \boldsymbol{C}_0

!     Transferring data to other subroutines
       statev(1)=H
       statev(2)=H

       if ((kflagS.eq.0).and.(npt.eq.8)) then
!     Computing gradient of damaged hydrostatic stress
       Var(1,1:8)=SH(locnum,1:8)
       do Nnpt=1,9
        call kGardVar(locnum,Nnpt,Var,GVar)
        do i=1,2
         GSH(locnum,Nnpt,i)=GVar(i)
        end do
       end do
      end if


       if (kflagS.eq.1) statev(2)=psit ! Staggered scheme


!!!!!! Part 2 including hydrogen concentration DOF
      elseif (cmname.eq.'MATERIAL-2') then

       stress=0.d0
       ddsdde=0.d0

      end if

      return
      end

c*****************************************************************

      subroutine umatht(u,dudt,dudg,flux,dfdt,dfdg,
     1statev,temp,dtemp,dtemdx,time,dtime,predef,dpred,
     2cmname,ntgrd,nstatv,props,nprops,coords,pnewdt,
     3noel,npt,layer,kspt,kstep,kinc)

      use ktransfer
      include 'aba_param.inc'

      character*80 cmname
      CHARACTER*80 CPNAME
      dimension dudg(ntgrd),flux(ntgrd),dfdt(ntgrd),
     1dfdg(ntgrd,ntgrd),statev(nstatv),dtemdx(ntgrd),
     2time(2),predef(1),dpred(1),props(nprops),coords(3)

      dimension GSHT(ntgrd)

!!!!!! Part 1 including mechanical and phase field DOFs
      if (CMNAME.eq.'MATERIAL-1') then

      CALL GETPARTINFO(noel, 1, CPNAME, locnum, JRCD)

      xl=props(1) ! Characteristic length scale, # Math: \ell
      Gc0=props(2) ! Fracture energy, # Math: G_c

       phi=temp+dtemp ! Phase field variable, # Math: \phi_{t+\delta t}=\phi_t+\delta \phi
       H=statev(2) ! Fracture driving force, # Math: \mathcal{H}
       cH=transcH(locnum,npt) ! Concentration of hydrogen, # Math: c_{\mathrm{H}}

       Vh=2000.d0 ! Molar volume of hydrogen, # Math: V_{\mathrm{H}}
       T=300.d0 ! Temperature, # Math: T
       R=8314.5d0 ! Gas constant, # Math: R

!     hydrogen contribution  
      Theta=cH*5.5d-05/(cH*5.5d-05+exp(-3.d7/(R*T))) ! Surface coverage, # Math: \theta=\frac{c_{\mathrm{H}}}{c_{\mathrm{H}}+\exp \left(\frac{-\Delta g_b^0}{R T_\mathrm{k}}\right)}
      Gc=(1.d0-0.89d0*Theta)*Gc0 ! Degraded fracture energy, # Math: G_c = (1 - \chi_{\mathrm{H}} \theta) G_{c0}

       U=U+(phi/xl**2-2.d0*(1.d0-phi)*H/(Gc*xl))*DTIME ! Internal energy, # Math: U_{t+\delta t} = U_t + \left(\frac{\phi}{\ell^2} - \frac{2(1-\phi)H}{G_c\ell}\right) \delta t
       DUDT=(1.d0/xl**2+2.d0*H/(Gc*xl))*DTIME ! Derivative of internal energy w.r.t. time, # Math: \frac{\partial U}{\partial \phi} = \frac{1}{\ell^2} + \frac{2H}{G_c \ell}
       DUDG=0.d0 ! Derivative of internal energy w.r.t. gradient of phase field, # Math: \frac{\partial U}{\partial \phi} = 0
       DFDT=0.d0 ! Derivative of flux w.r.t. phase field, # Math: \frac{\partial \mathbf{f}}{\partial \phi} = 0
       do i=1,NTGRD
        DFDG(i,i)=-1.d0 ! Derivative of flux w.r.t. gradient of phase field, # Math: \frac{\partial \mathbf{f}}{\partial \nabla \phi} = -\boldsymbol{I}
       end do
       FLUX=matmul(DFDG,DTEMDX) ! Flux, # Math: \mathbf{f} = - \nabla \phi

!!!!!! Part 2 including hydrogen concentration DOF
      elseif (CMNAME.eq.'MATERIAL-2') then

       CALL GETPARTINFO(noel, 1, CPNAME, locnum, JRCD) ! Get local element number

       D=props(1) ! Diffusion coefficient, # Math: D_{\mathrm{H}}

!     Updating hydrogen concentration
       cH=temp+dtemp ! Hydrogen concentration, # Math: c_{\mathrm{H}_{t+\delta t}}=c_{\mathrm{H}_t}+\delta c_{\mathrm{H}}

       Vh=2000.d0 ! Molar volume of hydrogen, # Math: V_{\mathrm{H}}
       T=300.d0 ! Temperature, # Math: T
       R=8314.5d0 ! Gas constant, # Math: R

!     Gradient of damaged hydrostatic stress
      do i=1,ntgrd
       GSHT(i)=GSH(locnum,npt,i) 
      end do

      u=u+dtemp ! Update internal energy, # Math: U_{t+\delta c_{\mathrm{H}}} = U_t + \delta c_{\mathrm{H}}
      dudt=1.d0 ! Derivative of internal energy w.r.t. time, # Math: \frac{\partial U}{\partial t} = 1
      dudg=0.0 ! Derivative of internal energy w.r.t. gradient of phase field, # Math: \frac{\partial U}{\partial \nabla \phi} = 0
      dfdt=(D*Vh*GSHT)/(R*T)
      do i=1,ntgrd
       dfdg(i,i)=-D ! Derivative of flux w.r.t. gradient of phase field, # Math: \frac{\partial \mathbf{f}}{\partial \nabla \phi} = -D \boldsymbol{I}
      end do
      flux=-D*dtemdx+cH*(D*Vh*GSHT)/(R*T) ! Flux, # Math: \mathbf{f} = -D \nabla c_{\mathrm{H}} + \frac{c_{\mathrm{H}} D V_{\mathrm{H}} \nabla \sigma_h}{R T}      

!     Transferring hydrogen concentration to UMAT subroutine
      transcH(locnum,npt)=cH

      end if

      return
      end

c*****************************************************************
!     Gradient of a variable
      subroutine kGardVar(locnum,npt,Var,GVar)

      use ktransfer
      include 'aba_param.inc'

      dimension Var(1,8), GVar(2)

      GVar=0.d0

      do i=1,8
       GVar(1)=GVar(1)+TransBx(locnum,npt,i)*Var(1,i)
       GVar(2)=GVar(2)+TransBy(locnum,npt,i)*Var(1,i)
      end do

      return
      end subroutine kGardVar

c*****************************************************************

      subroutine kshapefcn(npt,dNdz)
c
      include 'aba_param.inc'
c
      parameter (gaussCoord=sqrt(3.d0/5.d0))
      dimension dN(1,8),dNdz(2,8),coord24(2,9)

      data  coord24 /-1.d0, -1.d0, 
     20.d0, -1.d0, 
     31.d0, -1.d0, 
     4-1.d0, 0.d0, 
     50.d0, 0.d0, 
     61.d0, 0.d0, 
     7-1.d0, 1.d0, 
     80.d0, 1.d0, 
     91.d0, 1.d0/

!     2D 8-nodes
!     determine (g,h,r)
      g=coord24(1,npt)*gaussCoord
      h=coord24(2,npt)*gaussCoord

! !     shape functions 
!       dN(1,1)=-0.25d0*(1.d0-g)*(1.d0-h)*(1.d0+g+h) ! Shape function 1, # Math: N_1 = -\frac{1}{4}(1-g)(1-h)(1+g+h)
!       dN(1,2)=-0.25d0*(1.d0+g)*(1.d0-h)*(1.d0-g+h) ! Shape function 2, # Math: N_2 = -\frac{1}{4}(1+g)(1-h)(1-g+h)
!       dN(1,3)=-0.25d0*(1.d0+g)*(1.d0+h)*(1.d0-g-h) ! Shape function 3, # Math: N_3 = -\frac{1}{4}(1+g)(1+h)(1-g-h)
!       dN(1,4)=-0.25d0*(1.d0-g)*(1.d0+h)*(1.d0+g-h) ! Shape function 4, # Math: N_4 = -\frac{1}{4}(1-g)(1+h)(1+g-h)
!       dN(1,5)=0.5d0*(1.d0-g*g)*(1.d0-h) ! Shape function 5, # Math: N_5 = \frac{1}{2}(1-g^2)(1-h)
!       dN(1,6)=0.5d0*(1.d0+g)*(1.d0-h*h) ! Shape function 6, # Math: N_6 = \frac{1}{2}(1+g)(1-h^2)
!       dN(1,7)=0.5d0*(1.d0-g*g)*(1.d0+h) ! Shape function 7, # Math: N_7 = \frac{1}{2}(1-g^2)(1+h)
!       dN(1,8)=0.5d0*(1.d0-g)*(1.d0-h*h) ! Shape function 8, # Math: N_8 = \frac{1}{2}(1-g)(1-h^2)

!     derivative d(Ni)/d(g)
      dNdz(1,1)=0.25d0*(1.d0-h)*(2.d0*g+h) ! derivative d(Ni)/d(g), # Math: \frac{\partial N_1}{\partial g} = \frac{1}{4}(1-h)(2g+h)
      dNdz(1,2)=0.25d0*(1.d0-h)*(2.d0*g-h) ! derivative d(Ni)/d(g), # Math: \frac{\partial N_2}{\partial g} = \frac{1}{4}(1-h)(2g-h)
      dNdz(1,3)=0.25d0*(1.d0+h)*(2.d0*g+h) ! derivative d(Ni)/d(g), # Math: \frac{\partial N_3}{\partial g} = \frac{1}{4}(1+h)(2g+h)
      dNdz(1,4)=0.25d0*(1.d0+h)*(2.d0*g-h) ! derivative d(Ni)/d(g), # Math: \frac{\partial N_4}{\partial g} = \frac{1}{4}(1+h)(2g-h)
      dNdz(1,5)=-g*(1.d0-h) ! derivative d(Ni)/d(g), # Math: \frac{\partial N_5}{\partial g} = -g(1-h)
      dNdz(1,6)=0.5d0*(1.d0-h*h) ! derivative d(Ni)/d(g), # Math: \frac{\partial N_6}{\partial g} = \frac{1}{2}(1-h^2)
      dNdz(1,7)=-g*(1.d0+h) ! derivative d(Ni)/d(g), # Math: \frac{\partial N_7}{\partial g} = -g(1+h)
      dNdz(1,8)=-0.5d0*(1.d0-h*h) ! derivative d(Ni)/d(g), # Math: \frac{\partial N_8}{\partial g} = -\frac{1}{2}(1-h^2)

!     derivative d(Ni)/d(h)
      dNdz(2,1)=0.25d0*(1.d0-g)*(g+2.d0*h) ! derivative d(Ni)/d(h), # Math: \frac{\partial N_1}{\partial h} = \frac{1}{4}(1-g)(g+2h)
      dNdz(2,2)=0.25d0*(1.d0+g)*(2.d0*h-g) ! derivative d(Ni)/d(h), # Math: \frac{\partial N_2}{\partial h} = \frac{1}{4}(1+g)(2h-g)
      dNdz(2,3)=0.25d0*(1.d0+g)*(2.d0*h+g) ! derivative d(Ni)/d(h), # Math: \frac{\partial N_3}{\partial h} = \frac{1}{4}(1+g)(2h+g)
      dNdz(2,4)=0.25d0*(1.d0-g)*(2.d0*h-g) ! derivative d(Ni)/d(h), # Math: \frac{\partial N_4}{\partial h} = \frac{1}{4}(1-g)(2h-g)
      dNdz(2,5)=-0.5d0*(1.d0-g*g) ! derivative d(Ni)/d(h), # Math: \frac{\partial N_5}{\partial h} = -\frac{1}{2}(1-g^2)
      dNdz(2,6)=-(1.d0+g)*h ! derivative d(Ni)/d(h), # Math: \frac{\partial N_6}{\partial h} = -(1+g)h
      dNdz(2,7)=0.5d0*(1.d0-g*g) ! derivative d(Ni)/d(h), # Math: \frac{\partial N_7}{\partial h} = \frac{1}{2}(1-g^2)
      dNdz(2,8)=-(1.d0-g)*h ! derivative d(Ni)/d(h), # Math: \frac{\partial N_8}{\partial h} = -(1-g)h

      return
      end
c*****************************************************************
      subroutine kjacobian(noel,npt,dNg_1)
!     Notation: djac - Jac determinant; xjaci - inverse of Jac matrix 
!     dNdx - shape functions derivatives w.r.t. global coordinates
      use ktransfer
      include 'aba_param.inc'

      dimension xjac(2,2),xjaci(2,2),coords(2,8),
     1dNdz(2,8),dNdx(2,8),dNg_1(8,8),BSH(2,8)

      xjac=0.d0 ! Initialize jacobian matrix
      nnode=8 ! Number of nodes in the element
      ndim=2 ! Number of dimensions

      call kshapefcn(npt,dNdz) ! Get shape functions and their derivatives

!     Nodal coordinates of the element
      do i=1,nnode
       coords(1,i)=CoordNode(ElemConnec(noel,i),1) ! X-coordinate of node i
       coords(2,i)=CoordNode(ElemConnec(noel,i),2) ! Y-coordinate of node i
      end do

!     Compute Jacobian matrix
      do inod=1,nnode 
       do idim=1,ndim 
        do jdim=1,ndim 
         xjac(jdim,idim)=xjac(jdim,idim)+
     1   dNdz(jdim,inod)*coords(idim,inod) 
        end do
       end do
      end do

      djac=xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1) ! Jacobian determinant, # Math: J = \det(\mathbf{J}) = J_{11} J_{22} - J_{12} J_{21}
      if (djac.gt.0.d0) then ! Jacobian is positive - o.k.
       xjaci(1,1)=xjac(2,2)/djac ! Inverse of jacobian matrix, # Math: J^{-1}_{11} = \frac{J_{22}}{J}
       xjaci(2,2)=xjac(1,1)/djac ! Inverse of jacobian matrix, # Math: J^{-1}_{22} = \frac{J_{11}}{J}
       xjaci(1,2)=-xjac(1,2)/djac ! Inverse of jacobian matrix, # Math: J^{-1}_{12} = -\frac{J_{12}}{J}
       xjaci(2,1)=-xjac(2,1)/djac ! Inverse of jacobian matrix, # Math: J^{-1}_{21} = -\frac{J_{21}}{J}
      else ! Negative or zero jacobian
       write(7,*)'WARNING: element',jelem,'has negative Jacobian'
      endif

      dNdx=matmul(xjaci,dNdz) ! Shape function derivatives w.r.t. global coordinates, # Math: \frac{\partial N_i}{\partial x_j} = J^{-1}_{ij} \frac{\partial N_i}{\partial g_j}

      BSH=matmul(dNdx,dNg_1) ! 

      do i=1,8
       TransBx(noel,npt,i)=BSH(1,i)
       TransBy(noel,npt,i)=BSH(2,i)
      end do

      return
      end
c*****************************************************************
!     Calculating Inverse of shape function matrix for integration points 1-8
      subroutine kInvNg(dNg_1)

      include 'aba_param.inc'

      dimension dNg_1(8,8)

      dNg_1(1,1) = 5.0d0 * sqrt(15.0d0) / 9.0d0 
      dNg_1(1,2) = 5.0d0 / 3.0d0 - 2.0d0 * sqrt(15.0d0) / 3.0d0
      dNg_1(1,3) = 5.0d0 * sqrt(15.0d0) / 18.0d0 - 5.0d0 / 6.0d0
      dNg_1(1,4) = 5.0d0 / 3.0d0 - 2.0d0 * sqrt(15.0d0) / 3.0d0
      dNg_1(1,5) = 10.0d0 * sqrt(15.0d0) / 9.0d0 - 4.0d0
      dNg_1(1,6) = 5.0d0 / 3.0d0 - 4.0d0 * sqrt(15.0d0) / 9.0d0
      dNg_1(1,7) = 5.0d0 * sqrt(15.0d0) / 18.0d0 - 5.0d0 / 6.0d0
      dNg_1(1,8) = 5.0d0 / 3.0d0 - 4.0d0 * sqrt(15.0d0) / 9.0d0

      dNg_1(2,1) = 0.0d0
      dNg_1(2,2) = -sqrt(15.0d0) / 9.0d0
      dNg_1(2,3) = 5.0d0 * sqrt(15.0d0) / 18.0d0 + 5.0d0 / 6.0d0
      dNg_1(2,4) = sqrt(15.0d0) / 9.0d0
      dNg_1(2,5) = -2.0d0 / 3.0d0
      dNg_1(2,6) = -sqrt(15.0d0) / 9.0d0
      dNg_1(2,7) = 5.0d0 / 6.0d0 - 5.0d0 * sqrt(15.0d0) / 18.0d0
      dNg_1(2,8) = sqrt(15.0d0) / 9.0d0

      dNg_1(3,1) = -5.0d0 * sqrt(15.0d0) / 9.0d0
      dNg_1(3,2) = 2.0d0 * sqrt(15.0d0) / 3.0d0 + 5.0d0 / 3.0d0
      dNg_1(3,3) = -5.0d0 * sqrt(15.0d0) / 18.0d0 - 5.0d0 / 6.0d0
      dNg_1(3,4) = 2.0d0 * sqrt(15.0d0) / 3.0d0 + 5.0d0 / 3.0d0
      dNg_1(3,5) = -10.0d0 * sqrt(15.0d0) / 9.0d0 - 4.0d0
      dNg_1(3,6) = 4.0d0 * sqrt(15.0d0) / 9.0d0 + 5.0d0 / 3.0d0
      dNg_1(3,7) = -5.0d0 * sqrt(15.0d0) / 18.0d0 - 5.0d0 / 6.0d0
      dNg_1(3,8) = 4.0d0 * sqrt(15.0d0) / 9.0d0 + 5.0d0 / 3.0d0

      dNg_1(4,1) = 0.0d0
      dNg_1(4,2) = sqrt(15.0d0) / 9.0d0
      dNg_1(4,3) = 5.0d0 / 6.0d0 - 5.0d0 * sqrt(15.0d0) / 18.0d0
      dNg_1(4,4) = -sqrt(15.0d0) / 9.0d0
      dNg_1(4,5) = -2.0d0 / 3.0d0
      dNg_1(4,6) = sqrt(15.0d0) / 9.0d0
      dNg_1(4,7) = 5.0d0 * sqrt(15.0d0) / 18.0d0 + 5.0d0 / 6.0d0
      dNg_1(4,8) = -sqrt(15.0d0) / 9.0d0

      dNg_1(5,1) = 0.0d0
      dNg_1(5,2) = sqrt(15.0d0) / 6.0d0 + 5.0d0 / 6.0d0
      dNg_1(5,3) = 0.0d0
      dNg_1(5,4) = 0.0d0
      dNg_1(5,5) = -2.0d0 / 3.0d0
      dNg_1(5,6) = 0.0d0
      dNg_1(5,7) = 0.0d0
      dNg_1(5,8) = 5.0d0 / 6.0d0 - sqrt(15.0d0) / 6.0d0

      dNg_1(6,1) = 0.0d0
      dNg_1(6,2) = 0.0d0
      dNg_1(6,3) = 0.0d0
      dNg_1(6,4) = 5.0d0 / 6.0d0 - sqrt(15.0d0) / 6.0d0
      dNg_1(6,5) = -2.0d0 / 3.0d0
      dNg_1(6,6) = sqrt(15.0d0) / 6.0d0 + 5.0d0 / 6.0d0
      dNg_1(6,7) = 0.0d0
      dNg_1(6,8) = 0.0d0

      dNg_1(7,1) = 0.0d0
      dNg_1(7,2) = 5.0d0 / 6.0d0 - sqrt(15.0d0) / 6.0d0
      dNg_1(7,3) = 0.0d0
      dNg_1(7,4) = 0.0d0
      dNg_1(7,5) = -2.0d0 / 3.0d0
      dNg_1(7,6) = 0.0d0
      dNg_1(7,7) = 0.0d0
      dNg_1(7,8) = sqrt(15.0d0) / 6.0d0 + 5.0d0 / 6.0d0

      dNg_1(8,1) = 0.0d0
      dNg_1(8,2) = 0.0d0
      dNg_1(8,3) = 0.0d0
      dNg_1(8,4) = sqrt(15.0d0) / 6.0d0 + 5.0d0 / 6.0d0
      dNg_1(8,5) = -2.0d0 / 3.0d0
      dNg_1(8,6) = 5.0d0 / 6.0d0 - sqrt(15.0d0) / 6.0d0
      dNg_1(8,7) = 0.0d0
      dNg_1(8,8) = 0.0d0

      end subroutine kInvNg

c*****************************************************************
!     Extracting nodal coordinates and element-to-node connectivity

      subroutine KNodeElemInfo(JOBNAME)

      use ktransfer
      include 'aba_param.inc'

      character*262 input
      character*256 JOBNAME, line
      dimension n(8)

!      Input file name
      ! If you don't have permision to access scratch file, you sould put the address of directory including the input file in the following line.
      input=trim(JOBNAME)//'_f.inp'
      open(newunit=iunit, file=input, status='old',action='read', iostat=ierr)

      ! Read each line until "*Node" is found
      do
       read(iunit, '(A)') line
       if (line(1:5) == "*Node") then
        exit ! Exit the loop when "*Node" is found
       endif
      end do

      ! Read nodal coordinates until "*Element" is found
      numnode=0
      do
       read(iunit, '(A)') line
       if (line(1:8) == "*Element") then
        exit ! Exit the loop when "*Element" is found
       endif
       numnode=numnode+1
       ! Extract nodal coordinates from the line
       read(line, *) nnode, x, y

       CoordNode(nnode,1)=x
       CoordNode(nnode,2)=y

      end do

      ! Close the input file
      close(iunit)

      ! Open the Abaqus input file
      open(newunit=iunit, file=input, status='old',action='read', iostat=ierr)

      ! Read each line until "*Node" is found
      do
       read(iunit, '(A)') line
       if (line(1:8) == "*Element") then
        exit ! Exit the loop when "*Node" is found
       endif
      end do

      ! Read nodal coordinates until "*Element" is found
      numElem=0
      do
       read(iunit, '(A)') line
       if (line(1:5) == "*Nset") then
        exit ! Exit the loop when "*Element" is found
       endif
       numElem=numElem+1
       ! Extract nodal coordinates from the line
       read(line, *) nElem, n(1), n(2), n(3), n(4), n(5), n(6), n(7), n(8) 

       do i=1,8
        ElemConnec(nElem,i)=n(i)
       end do

      end do

!     Close the input file
      close(iunit)

      return
      end subroutine KNodeElemInfo
c*****************************************************************




