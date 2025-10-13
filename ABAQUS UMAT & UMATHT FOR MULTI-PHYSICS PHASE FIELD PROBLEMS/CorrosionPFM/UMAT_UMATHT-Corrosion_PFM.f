
!----------------------------------------------------------------------
! UMAT and UMATHT Subroutines for Phase Field Corrosion Modelling
!----------------------------------------------------------------------
!
! This code implements a simple phase field approach for modelling
! corrosion processes. It includes user subroutines UMAT and UMATHT
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
! "A phase field model for hydraulic fracture: Drucker–Prager driving force and a hybrid coupling strategy",
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
      module ktansfer

      real*8 :: TranscL(700000,4) = 0.d0
      real*8 :: TransDphiDx(700000,4,3) = 0.d0
      real*8 :: TransDhphi(700000,4) = 0.d0
      real*8 xkap, cse, cle, D

      save
      end module

c*****************************************************************

      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,predef,dpred,
     2cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      use ktansfer
      include 'aba_param.inc'

      character*80 cmname
      CHARACTER*80 CPNAME
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),time(2),
     2predef(1),dpred(1),props(nprops),coords(3),drot(3,3),dfgrd0(3,3),
     3dfgrd1(3,3),jstep(4)

      dimension eelas(ntens),eplas(ntens),flow(ntens),olds(ntens),
     2oldpl(ntens)

      parameter(toler=1.d-6,newton=20)

!     Mechanical and Phase field part
      if (cmname.eq.'MATERIAL-1') then

       CALL GETPARTINFO(noel, 1, CPNAME, locnum, JRCD)

       ! Elastic strain, # Math: \boldsymbol{\varepsilon}_{e}
       eelas(1:ntens)=statev((2*ntens+1):3*ntens)
       ! Plastic strain, # Math: \boldsymbol{\varepsilon}_{p}
       eplas(1:ntens)=statev((3*ntens+1):4*ntens)
       ! Equivalent plastic strain, # Math: \varepsilon^{p}
       eqplas=statev(1+4*ntens)
       ! Increment of equivalent plastic strain, # Math: \delta \varepsilon^{p}
       deqpl=statev(5+4*ntens)
       olds=statev(1:ntens)
       oldpl=eplas

!     Initialization
       ddsdde=0.d0
       E=props(1) ! Young's modulus, # Math: E
       xnu=props(2) ! Poisson's ratio, # Math: \nu
       Sy=props(3) ! Yield stress, # Math: \sigma_{y}
       xn=props(4) ! Strain hardening exponent, # Math: n
       xkap=1.d-7  ! well-conditioning parameter

!     Updating and transferring variables at the beginning of increment
       phi=temp+dtemp ! Phase field, # Math: \phi
       hphi=-2.d0*phi**3+3.d0*phi**2 ! Phase field interpolation function, # Math: g(\phi)

!     Build elastic stiffness matrix
       eg=E/(1.d0+xnu)/2.d0 ! Shear modulus, # Math: \mu= \frac{E}{2(1+\nu)}
       elam=(E/(1.d0-2.d0*xnu)-2.d0*eg)/3.d0 ! Lame's first parameter, # Math: \lambda= \frac{E \nu}{(1+\nu)(1-2\nu)}

       do i=1,3
        do j=1,3
         ddsdde(j,i)=elam
        end do
        ddsdde(i,i)=2.d0*eg+elam
       end do
       do i=4,ntens
        ddsdde(i,i)=eg
       end do


!     Calculate predictor stress and elastic strain
       stress=olds+matmul(ddsdde,dstran) ! Predictor stress, # Math: \boldsymbol{\sigma}=\boldsymbol{\sigma}^{old}+\boldsymbol{C} \delta \boldsymbol{\varepsilon}_{e}
       eelas=eelas+dstran ! Elastic strain, # Math: \boldsymbol{\varepsilon}_{e}=\boldsymbol{\varepsilon}_{e}^{old}+\delta \boldsymbol{\varepsilon}_{e}

!     Calculate equivalent von Mises stress
       Smises=(stress(1)-stress(2))**2+(stress(2)-stress(3))**2 
     1 +(stress(3)-stress(1))**2 
       do i=4,ntens
        Smises=Smises+6.d0*stress(i)**2
       end do
       Smises=sqrt(Smises/2.d0) ! Mises stress, # Math: \sigma_{m}=\sqrt{\frac{1}{2}((\sigma_{1}-\sigma_{2})^{2}+(\sigma_{2}-\sigma_{3})^{2}+(\sigma_{3}-\sigma_{1})^{2}+3(\sigma_{4}^{2}+\sigma_{5}^{2}+\sigma_{6}^{2}))}

!     Get yield stress from the specified hardening curve
       Sf=Sy*(1.d0+E*eqplas/Sy)**xn ! Flow stress, # Math: \sigma_{f}=\sigma_{y}(1+\frac{E}{\sigma_{y}} \varepsilon^{p})^{n}

!     Determine if active yielding
       if (Smises.gt.(1.d0+toler)*Sf) then

!     Calculate the flow direction
        Sh=(stress(1)+stress(2)+stress(3))/3.d0 ! Hydrostatic stress, # Math: \sigma_{h}=\frac{1}{3}(\sigma_{1}+\sigma_{2}+\sigma_{3})
        flow(1:3)=(stress(1:3)-Sh)/Smises
        flow(4:ntens)=stress(4:ntens)/Smises

!     Solve for Smises and deqpl using Newton's method
        deqpl=0.d0
        Et=E*xn*(1.d0+E*eqplas/Sy)**(xn-1) ! Derivative of flow stress, # Math: \frac{d\sigma_{f}}{d\varepsilon^{p}}=\frac{E n}{\sigma_{y}}(1+\frac{E}{\sigma_{y}} \varepsilon^{p})^{n-1}
        do kewton=1,newton
         rhs=Smises-(3.d0*eg)*deqpl-Sf ! Residual, # Math: \sigma_{m}-(3\mu \delta \varepsilon^{p}+\sigma_{f})
         deqpl=deqpl+rhs/((3.d0*eg)+Et) ! Update increment of equivalent plastic strain, # Math: \delta \varepsilon^{p}=\frac{\sigma_{m}-\sigma_{f}}{3\mu+E n (1+\frac{E}{\sigma_{y}} \varepsilon^{p})^{n-1}}
!        if(deqpl.lt.0.d0) deqpl=-deqpl
         Sf=Sy*(1.d0+E*(eqplas+deqpl)/Sy)**xn ! Update flow stress, # Math: \sigma_{f}=\sigma_{y}(1+\frac{E}{\sigma_{y}} (\varepsilon^{p}+\delta \varepsilon^{p}))^{n}
         Et=E*xn*(1.d0+E*(eqplas+deqpl)/Sy)**(xn-1) ! Update derivative of flow stress, # Math: \frac{d\sigma_{f}}{d\varepsilon^{p}}=\frac{E n}{\sigma_{y}}(1+\frac{E}{\sigma_{y}} (\varepsilon^{p}+\delta \varepsilon^{p}))^{n-1}
         if(abs(rhs).lt.toler*Sy) exit
        end do

        if (kewton.eq.newton) write(7,*)'WARNING: plasticity loop failed'

! Update stresses and strains
        stress(1:3)=flow(1:3)*Sf+Sh 
        eplas(1:3)=eplas(1:3)+3.d0/2.d0*flow(1:3)*deqpl
        eelas(1:3)=eelas(1:3)-3.d0/2.d0*flow(1:3)*deqpl
        stress(4:ntens)=flow(4:ntens)*Sf
        eplas(4:ntens)=eplas(4:ntens)+3.d0*flow(4:ntens)*deqpl
        eelas(4:ntens)=eelas(4:ntens)-3.d0*flow(4:ntens)*deqpl 
        eqplas=eqplas+deqpl ! Equivalent plastic strain, # Math: \varepsilon^{p}=\varepsilon^{p}+\delta \varepsilon^{p}

!    Calculate the plastic strain energy density
        do i=1,ntens
         spd=spd+(stress(i)+olds(i))*(eplas(i)-oldpl(i))/2.d0
        end do

!     Formulate the Jacobian (material tangent)   
        effg=eg*Sf/Smises
        efflam=(E/(1.d0-2.d0*xnu)-2.d0*effg)/3.d0
        effhrd=3.d0*eg*Et/(3.d0*eg+Et)-3.d0*effg
        do i=1,3
         do j=1,3
          ddsdde(j,i)=efflam
         enddo
         ddsdde(i,i)=2.d0*effg+efflam
        end do
        do i=4,ntens
         ddsdde(i,i)=effg
        end do

        do i=1,ntens
         do j=1,ntens
          ddsdde(j,i)=ddsdde(j,i)+effhrd*flow(j)*flow(i) 
         end do
        end do
       endif


!     Updating and transferring variables at the end of increment
       Sh=(stress(1)+stress(2)+stress(3))/3.d0
       statev(1:ntens)=stress(1:ntens)
       statev((2*ntens+1):3*ntens)=eelas
       statev((3*ntens+1):4*ntens)=eplas
       statev(4*ntens+1)=eqplas
       statev(4*ntens+4)=Sh
       statev(4*ntens+5)=deqpl


       stress=(hphi+xkap)*stress
       ddsdde=(hphi+xkap)*ddsdde

!     Mass diffusion part
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

      use ktansfer
      include 'aba_param.inc'

      character*80 cmname
      CHARACTER*80 CPNAME
      dimension dudg(ntgrd),flux(ntgrd),dfdt(ntgrd),
     1dfdg(ntgrd,ntgrd),statev(nstatv),dtemdx(ntgrd),
     2time(2),predef(1),dpred(1),props(nprops),coords(3)

      dimension DphiDX(ntgrd), FF(ntgrd)

!     Mechanical and Phase field part
      if (cmname.eq."MATERIAL-1") then

       CALL GETPARTINFO(noel, 1, CPNAME, locnum, JRCD)


       !     Reading parameters
       E=props(1) ! Young's modulus, # Math: E
       xnu=props(2) ! Poisson's ratio, # Math: \nu
       Sy=props(3) ! Yield stress, # Math: \sigma_y
       xn=props(4) ! Strain hardening exponent, # Math: n
       D=props(5) ! Diffusion coefficient, # Math: D
       xL0=props(6) ! Initial interface kinetics coefficient, # Math: L_0
       a_phi=props(7) ! Gradient energy coefficient, # Math: \kappa
       wh=props(8) ! Height of the double well potential, # Math: \omega
       AA=props(9) ! Free energy density curvature, # Math: A
       xk=props(10) ! Well-conditioning parameter
       ef=props(11) 
       t0=props(12) ! Time constant for the decay of interface kinetics
       cse=1.d0 ! Normalized equilibrium concentration of solid
       cle=0.036d0 ! Normalized equilibrium concentration of liquid
      !  xkap=1.d-7 ! well-conditioning parameter
       T=300.d0 ! Temperature (K), # Math: T
       R=8314.d0 ! Gas constant, # Math: R
       xkap=1.d-7  ! Well-conditioning parameter


!     Updating and transferring variables at the beginning of increment
       phi=temp+dtemp
       ntens=2*ntgrd ! Number of stress components
       pls=statev(4*ntens+1) ! Plastic equivalent strain
       Sh0=statev(4*ntens+4) ! Undamaged hydrostatic stress
       ti=statev(4*ntens+7) ! Accumulated time
       ei=statev(4*ntens+8) ! Corrosion density
       cL=TranscL(locnum,npt) ! Concentration


!     Defining suitable functions       
       hphi=-2.d0*phi**3+3.d0*phi**2 ! Phase field interpolation function, # Math: g(\phi)=-2\phi^{3}+3\phi^{2}
       dhphi=-6.d0*phi**2+6.d0*phi ! Derivative of phase field interpolation function, # Math: \frac{\partial g(\phi)}{\partial\phi}=-6\phi^{2}+6\phi
       ddhphi=-12.d0*phi+6.d0 ! Second derivative of phase field interpolation function, # Math: \frac{\partial^{2} g(\phi)}{\partial\phi^{2}}=-12\phi+6
       gphi=(1.d0-phi)*(1.d0-phi)*phi**2 ! Double well potential, # Math: w(\phi)=(1-\phi)^{2}\phi^{2}
       dgphi=2.d0*phi+4.d0*phi**3-6.d0*phi**2 ! Derivative of double well potential, # Math: \frac{\partial w(\phi)}{\partial\phi}=2\phi+4\phi^{3}-6\phi^{2}
       ddgphi=2.d0+12.d0*phi**2-12.d0*phi ! Second derivative of double well potential, # Math: \frac{\partial^{2} w(\phi)}{\partial\phi^{2}}=2+12\phi^{2}-12\phi


!     Enhance the value of L
       Sh=(hphi+xkap)*Sh0 ! Degradation of hydrostatic stress, # Math: (g(\phi)+\kappa)\sigma_{h0}
       xkm=exp(Sh*7.12d3/(R*T))*(1.d0+pls/(props(3)/props(1))) 
       if ((time(2)+dtime).lt.3.01d0) then
        ti=0.d0
        ei=0.d0
       else
        ei=ei+statev(4*ntens+5)
        ti=ti+dtime
       endif

       if (ei.gt.ef) then
        ti=0.d0
        ei=0.d0
       endif

       if (ti.lt.t0) then
        xL=xL0*xkm
       else
        xL=xL0*xkm*exp(-xk*(ti-t0))
       endif

       dfc=-2.d0*AA*(cL-hphi*(cse-cle)-cle)*(cse-cle)*dhphi+wh*dgphi
       ddfc=wh*ddgphi-2.d0*AA*(cse-cle)*
     1 (ddhphi*(cL-hphi*(cse-cle)-cle)-dhphi*dhphi*(cse-cle))


       U=U+(dtemp/xL+dfc*DTIME)
       DUDT=(1.d0/xL+ddfc*DTIME)
       DUDG=0.d0
       DFDT=0.d0
       do i=1,NTGRD
        DFDG(i,i)=-a_phi
       end do
       FLUX=matmul(DFDG,DTEMDX)


!     Updating and transferring variables at the end of increment
       statev(4*ntens+7)=ti
       statev(4*ntens+8)=ei
       do i=1,ntgrd
        TransDphiDx(locnum,npt,i)=DTEMDX(i)
       end do
       TransDhphi(locnum,npt)=dhphi

!     Mass diffusion part
      elseif (cmname.eq."MATERIAL-2") then

       CALL GETPARTINFO(noel, 1, CPNAME, locnum, JRCD)



!     Updating and transferring variables at the beginning of increment
       cL=temp+dtemp ! Updating concentration
       do i=1,ntgrd ! Importing gradient of phase field
        DphiDX(i)=TransDphiDx(locnum,npt,i)
       end do
       dhphi=TransDhphi(locnum,npt)

       FF=DTEMDX-(cse-cle)*dhphi*DphiDX

       U=U+dtemp
       DUDT=1.d0
       DUDG=0.d0
       DFDT=0.d0
       do i=1,NTGRD
        DFDG(i,i)=-D
       end do
       FLUX=matmul(DFDG,FF)

!     Updating and transferring variables at the end of increment
       TranscL(locnum,npt)=cL ! Exporting concentration

      end if

      return
      end

c*****************************************************************
