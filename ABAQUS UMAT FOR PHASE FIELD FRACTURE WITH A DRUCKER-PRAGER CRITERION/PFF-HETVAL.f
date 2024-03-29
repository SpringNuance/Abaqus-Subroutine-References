! UMAT and HETVAL subroutines for implementing phase field fracture
! The code is distributed under a BSD license     

! If using this code for research or industrial purposes, please cite:
! Y. Navidtehrani, C. Betegon, E. Martinez-Paneda. A simple and robust
! Abaqus implementation of the phase field fracture method
! Applications in Engineering Science 6, 100050 (2021)
! https://doi.org/10.1016/j.apples.2021.100050
! &
! Navidtehrani, Y.; Betegón, C.; Martínez-Pañeda, E. A general framework
! for decomposing the phase field fracture driving force, particularised
! to a Drucker–Prager failure surface. Theoretical and Applied Fracture
! Mechanics 2022, 121, 103555. 
! https://doi.org/10.1016/j.tafmec.2022.103555

! Emilio Martinez-Paneda (e.martinez-paneda@imperial.ac.uk)
! Imperial College London

c***********************************************************************

      subroutine hetval(cmname,temp,time,dtime,statev,flux,predef,dpred)

      include 'aba_param.inc'

      character*80 cmname

      dimension temp(2),statev(*),predef(*),time(2),flux(2),dpred(*)

      real*8, parameter :: pi = 3.1415926535897932384626433832d0

      phi=temp(1)

      H=statev(2)
      xl=statev(3)
      Gc=statev(4)
      kflagC=statev(5)
      dg=statev(6)
      ddg=statev(7)

      if (kflagC.eq.0) then
       w=phi**2
       dw=2.d0*phi
       ddw=2.d0
       cw=0.5d0
      elseif (kflagC.eq.1) then
       w=phi
       dw=1.d0
       ddw=0.d0
       cw=2.d0/3.d0
      else
       w=2.d0*phi-phi**2
       dw=2.d0-2.d0*phi
       ddw=-2.d0
       cw=pi/4.d0
      end if

      flux(1)=-(dg*H*2.d0*cw/(xl*Gc)+dw/(2.d0*xl**2))
      flux(2)=-(ddg*H*2.d0*cw/(xl*Gc)+ddw/(2.d0*xl**2))

      return
      end
c***********************************************************************

      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,predef,dpred,
     2cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      include 'aba_param.inc'

      character*80 cmname
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),time(2),
     2predef(1),dpred(1),props(nprops),coords(3),drot(3,3),dfgrd0(3,3),
     3dfgrd1(3,3),jstep(4)

      dimension PS(3), AN(3,3)

      real*8, parameter :: pi = 3.1415926535897932384626433832d0

!     Initialization
      ddsdde=0.d0
      Hmin=0.d0
      E=props(1)! ! Young's modulus
      xnu=props(2) ! Poisson's ratio
      xl=props(3)! ! Phase field length scale
      Gc=props(4) ! Toughness
      xk=1.d-07
      kflagS=int(props(5)) ! Solution flag (0: monolithic, 1: staggered)
      kflagC=int(props(6)) ! Model (0:AT2, 1:AT1, 2:PF-CZM/linear, 3:PF-CZM/exp)
      kflagD=int(props(7)) ! Split (0: No split, 1: Amor, 2: Miehe, 3: Drucker-Prager based)
      kflagH=int(props(8)) ! Split linear momentum (1:Hybrid/no split, 2:Anisotropic)
      if (kflagC.eq.2.or.kflagC.eq.3) then
       ft=props(9) ! Tensile strength
       a=4.d0*E*Gc/(pi*xl*ft**2)
      end if
      if (kflagD.eq.3) BDP=props(10) ! B parameter of Drucker-Prager model

      phi=temp+dtemp
      psit=statev(1)

      eg=E/(1.d0+xnu)/2.d0
      elam=(E/(1.d0-2.d0*xnu)-2.d0*eg)/3.d0
      bk=E/(1.d0-2.d0*xnu)/3.d0

!     Degradation function
      call KDegFun(kflagC,phi,xk,a,g,dg,ddg)

!     Build stiffness matrix
      call KStiffMat(E,xnu,eg,elam,bk,BDP,ndi,nshr,ntens,stran,ddsdde,
     1kflagH,kflagD,g)

!     Update strains and stresses      
      stran=stran+dstran
      stress=matmul(ddsdde,stran)


!     Compute the driving force for fracture
      call KDriFor(E,xnu,eg,elam,bk,BDP,Gc,xl,ft,ndi,nshr,ntens,stran,
     1stress,kflagC,kflagD,psit,g,H,Hmin)

!     State variables
      statev(1)=H
      statev(2)=H
      statev(3)=xl
      statev(4)=Gc
      statev(5)=kflagC
      statev(6)=dg
      statev(7)=ddg

      if (kflagS.eq.1) statev(2)=max(psit,Hmin) ! Staggered scheme

      return
      end
c***********************************************************************

!     Degradation function
      subroutine KDegFun(kflagC,phi,xk,a,g,dg,ddg)

      include 'aba_param.inc'

      ! Variables returned from the subroutine:
      ! (g,dg,ddg)

      ! Variables to be provided to the subroutine:
      ! (kflagC,phi,xk,a)

      if (kflagC.eq.0.or.kflagC.eq.1) then ! AT2 & AT1 model
       g=(1.d0-phi)**2+xk
       dg=-2.d0*(1.d0-phi)
       ddg=2.d0
      else ! PF-CZM model
       if (kflagC.eq.2) then ! Linear PF-CZM
        b=-0.5d0
        c=0.d0
        d=2.d0
       elseif (kflagC.eq.3) then ! Exponential PF-CZM
        b=2.d0**(5.d0/3.d0)-3.d0
        c=0.d0
        d=2.5d0
       endif
       gNum=(1.d0-phi)**d
       dgNum=-d*(1.d0-phi)**(d-1.d0);
       ddgNum=d*(d-1.d0)*(1.d0-phi)**(d-2.d0)
       gDen=gNum+a*phi+a*b*phi**2.d0+a*b*c*phi**3.d0
       dgDen=dgNum+a+2.d0*a*b*phi+3.d0*a*b*c*phi**2.d0
       ddgDen=ddgNum+2.d0*a*b+6.d0*a*b*c*phi

       g=gNum/gDen
       dg=(dgNum*gDen-gNum*dgDen)/(gDen**2.d0)
       ddg=((ddgNum*gDen-gNum*ddgDen)*gDen-2.d0*
     1 (dgNum*gDen-gNum*dgDen)*dgDen)/(gDen**3.d0)

      end if

      end
c***********************************************************************

!     Build stiffness matrix
      subroutine KStiffMat(E,xnu,eg,elam,bk,BDP,ndi,nshr,ntens,stran,ddsdde,
     1kflagH,kflagD,g)

      include 'aba_param.inc'

      dimension stress(ntens),ddsdde(ntens,ntens),stran(ntens),PS(3),
     1AN(3,3)

      ! Variables returned from the subroutine:
      ! (ddsdde)

      ! Variables to be provided to the subroutine:
      ! (E,xnu,eg,elam,ndi,nshr,ntens,bk,BDP,stran,kflagH,kflagD,g)

      eg=E/(1.d0+xnu)/2.d0
      elam=(E/(1.d0-2.d0*xnu)-2.d0*eg)/3.d0
      bk=E/(1.d0-2.d0*xnu)/3.d0

      if ((kflagH.eq.1).or.(all(stran.eq.0.d0))) then ! Hybrid method
       if (ndi .eq. 3) then ! 3D or Plane strain
        do i=1,3
         do j=1,3
          ddsdde(j,i)=elam
         end do
         ddsdde(i,i)=eg*2.d0+elam
        end do
        do i=4,ntens
         ddsdde(i,i)=eg
        end do
       elseif (ndi.eq.2) then ! Plane stress
        ddsdde(1,1)=E/(1.d0-xnu**2)
        ddsdde(1,2)=E*xnu/(1.d0-xnu**2)
        ddsdde(2,1)=E*xnu/(1.d0-xnu**2)
        ddsdde(2,2)=E/(1.d0-xnu**2)
        ddsdde(3,3)=eg
       endif
       ddsdde=ddsdde*g
      elseif (kflagH .eq. 2) then ! Anisotropic method
       call kAnisEla(stress,stran,ddsdde,ndi,nshr,ntens,g,eg,
     1 elam,bk,BDP,kflagD)
      endif

      end
c***********************************************************************

!     Calculating 4th-order anisotropic elasticity matrix 
      subroutine kAnisEla(stress,stran,ddsdde,ndi,nshr,ntens,g,eg,
     1elam,bk,BDP,kflagD)

      include 'aba_param.inc'

      dimension stress(ntens),stran(ntens),ddsdde(ntens,ntens),
     1ddsdde0(6,6),D3ddsdde(6,6),ddPsidE(6,6),ddPsipdE(6,6),
     2ddPsindE(6,6),PV(3),PVp(3),PVn(3),dVpdV(6),dVndV(6),AN(3,3),
     3PrStran(6),PrStranP(6),PrStranN(6),R(6,6),PS(3),PSs(3)

      dimension dI1dV(6),ddI1dV(6,6),dI1pdV(6),ddI1pdV(6,6),dI1ndV(6),
     1ddI1ndV(6,6),dI2dV(6),ddI2dV(6,6),dI2pdV(6),ddI2pdV(6,6),
     2dI2ndV(6),ddI2ndV(6,6),dI3dV(6),ddI3dV(6,6),dI3pdV(6),
     3ddI3pdV(6,6),dI3ndV(6),ddI3ndV(6,6),dI1Sds(6,1),stress0(ntens),
     4dJ2dV(6),ddJ2dV(6,6)

      ! Variables returned from the subroutine:
      ! (ddsdde)

      ! Variables to be provided to the subroutine:
      ! (stress,stran,ndi,nshr,ntens,g,eg,elam,bk,BDP,kflagD)

!     Initiation
      ddPsidE=0.d0 ! Elastic stiffness tensor
      ddPsipdE=0.d0 ! Damaged part of Elastic stiffness tensor 
      ddPsindE=0.d0 ! Undamaged part of Elastic stiffness tensor
      ddsdde0=0.d0 ! Elastic stiffness tensor of Effective configuration
      D3ddsdde=0.d0 ! 3D Elastic stiffness tensor
      PrStran=0.d0 ! Principal strain

!     Elastic stiffness tensor of Effective configuration
      do i=1,3
       do j=1,3
        ddsdde0(j,i)=elam
       end do
       ddsdde0(i,i)=eg*2.d0+elam
      end do
      do i=4,6
       ddsdde0(i,i)=eg
      end do

      if (kflagD .eq. 1) then ! Amor et al.

       trE=sum(stran(1:3)) ! Trace of strain

       if (trE .lt. 0.d0) then
        ddPsindE(1:3,1:3)=bk
       elseif (trE .eq. 0.d0) then
        ddPsindE(1:3,1:3)=5.d-1*bk
       end if

!     Damaged elasticity tensor 
       ddPsidE=g*ddsdde0+(1.d0-g)*ddPsindE
       D3ddsdde=ddPsidE

      else if (kflagD .eq. 2) then ! Miehe et al.
      
       call SPRIND(stran,PV,AN,2,ndi,nshr)

       PrStran(1:3)=PV ! principal strain

       PrStranP=(PrStran+abs(PrStran))/2.d0 ! Positive part of principal strain

       call kInvariant(PrStranP,xI1p,dI1pdV,ddI1pdV,1,ndi,nshr)
       call kInvariant(PrStranP,xI2p,dI2pdV,ddI2pdV,2,ndi,nshr)

!     Variation of positive part on original part
       dVpdV=0.d0
       do i=1,6
        if (PrStran(i) .gt. 0.d0) then
         dVpdV(i)=1.d0
        elseif (PrStran(i) .eq. 0.d0) then
         dVpdV(i)=5.d-1
        endif
       end do

       if (xI1 .gt. 0.d0) then
        ddPsipdE(1:3,1:3)=elam
       elseif (xI1 .eq. 0.d0) then
        ddPsipdE(1:3,1:3)=5.d-1*elam
       end if

       do i=1,6
        do j=1,6
         ddPsipdE(i,j)=ddPsipdE(i,j)+2.d0*eg*(dI1pdV(i)*dVpdV(i)*
     1   dI1pdV(j)*dVpdV(j)-ddI2pdV(i,j)*dVpdV(i)*dVpdV(j))
        end do
       end do

!     Damaged elasticity tensor in principal directions
       ddPsidE=(g-1.d0)*ddPsipdE+ddsdde0

!     Damaged elasticity tensor in the original directions
       call kRotation4thOrder(D3ddsdde,ddPsidE,transpose(AN))

      else if (kflagD .eq. 3) then ! Drucker-Prager based split, Navidtehrani et al.

!     1st and 2nd Invariants of strain and their derivations
       call kInvariant(stran,xI1,dI1dV,ddI1dV,1,ndi,nshr)
       call kInvariant(stran,xI2,dI2dV,ddI2dV,2,ndi,nshr)

!     2nd Invariant of deviatoric part of strain
       xJ2=xI1**2/3.d0-xI2
       dJ2dV=(2.d0/3.d0)*xI1*dI1dV-dI2dV
       do i=1,6
        do j=1,6
         ddJ2dV(i,j)=(2.d0/3.d0)*(dI1dV(i)*dI1dV(j))-ddI2dV(i,j)
        end do
       end do

       if (-6.d0*BDP*sqrt(xJ2).lt.xI1) then
        ddPsindE=0.d0
       elseif (2.d0*eg*sqrt(xJ2).lt.(3.d0*BDP*bk*xI1)) then
        ddPsindE=ddsdde0
       else
        do i=1,6
         do j=1,6
          ddPsindE(i,j)=(bk*eg/(9.d0*BDP**2*bk+eg))*
     1    (dI1dV(i)+3.d0*BDP*dJ2dV(i)/sqrt(xJ2))*
     2    (dI1dV(j)+3.d0*BDP*dJ2dV(j)/sqrt(xJ2))+
     3    6.d0*BDP*(bk*eg/(18.d0*BDP**2*bk+2.d0*eg))*
     4    (xI1+6.d0*BDP*sqrt(xJ2))/sqrt(xJ2)*
     5    (ddJ2dV(i,j)-dJ2dV(i)*dJ2dV(j)/xJ2/2.d0)
         end do
        end do
       endif
       ddPsidE=g*ddsdde0+(1.d0-g)*ddPsindE
       D3ddsdde=ddPsidE
      end if


      ddsdde=D3ddsdde(1:ntens,1:ntens)


      end subroutine kAnisEla

c***********************************************************************

!     Compute the driving force for fracture
      subroutine KDriFor(E,xnu,eg,elam,bk,BDP,Gc,xl,ft,ndi,nshr,ntens,stran
     1,stress,kflagC,kflagD,psit,g,H,Hmin)

      include 'aba_param.inc'

      dimension stress(ntens),ddsdde(ntens,ntens),stran(ntens),
     1Edev(ntens),PS(3),AN(3,3),PSs(3)

      ! Variables returned from the subroutine:
      ! (H)

      ! Variables to be provided to the subroutine:
      ! (E,xnu,eg,elam,bk,BDP,Gc,xl,ft,ndi,nshr,ntens,stran,stress,kflagC,kflagD,psit,g,Hmin)

      if (kflagC.eq.0.or.kflagC.eq.1) then
       psip=0.d0

       if (kflagD.eq.1) then ! Amor et al.
        CALL SPRINC(stran,PS,2,ndi,nshr)
        Edev=PS
        trE=PS(1)+PS(2)+PS(3)
        Edev(1:3)=PS(1:3)-trE/3.d0
        EdevS=Edev(1)**2+Edev(2)**2+Edev(3)**2
        trEp=0.5d0*(trE+abs(trE))
        trEn=0.5d0*(trE-abs(trE))
        psip=0.5d0*bk*trEp**2+eg*EdevS
        psin=0.5d0*bk*trEn**2

       elseif (kflagD.eq.2) then ! Miehe et al.
        call SPRINC(stran,PS,2,ndi,nshr)
        trp1=(PS(1)+PS(2)+PS(3)+abs(PS(1)+PS(2)+PS(3)))/2.d0
        trn1=(PS(1)+PS(2)+PS(3)-abs(PS(1)+PS(2)+PS(3)))/2.d0
        trp2=0.d0
        trn2=0.d0
        do i=1,3
         trp2=trp2+(PS(i)+abs(PS(i)))**2.d0/4.d0
         trn2=trn2+(PS(i)-abs(PS(i)))**2.d0/4.d0
        end do
        psip=xnu*eg/(1d0-2d0*xnu)*trp1**2d0+eg*trp2
        psin=xnu*eg/(1d0-2d0*xnu)*trn1**2d0+eg*trn2

       else  if (kflagD.eq.3) then ! Drucker-Prager based split, Navidtehrani et al.

        CALL SPRINC(stran,PS,2,ndi,nshr)

        xI1=sum(PS(1:3)) ! 1st invariant of strain tensor
        xJ2=((PS(1)-PS(2))**2+(PS(2)-PS(3))**2+(PS(3)-PS(1))**2)/6.d0 ! 2nd invariant of deviatoric part of tensor

        if (-6.d0*BDP*sqrt(xJ2).lt.xI1) then
         psip=0.5d0*bk*xI1**2+2.d0*eg*xJ2
        elseif (2.d0*eg*sqrt(xJ2).lt.(3.d0*BDP*bk*xI1)) then
         psip=0.d0
        else
         psip=(-3.d0*BDP*bk*xI1+2.d0*eg*sqrt(xJ2))**2/
     1   (18.d0*BDP**2*bk+2.d0*eg)
        endif

       else ! no split
        psip=0.d0
        do i=1,ntens
         psip=psip+stress(i)*stran(i)*0.5d0/g
        end do
       endif

       if (kflagC.eq.1) Hmin=3.d0*Gc/(16.d0*xl)
       H=max(psit,psip,Hmin)
      else
       call SPRINC((stress/g),PS,1,ndi,nshr)
       psip=0.5d0*max(PS(1),PS(2),PS(3))**2/E
       Hmin=0.5d0*ft**2/E
       H=max(psit,psip,Hmin)
      end if

      end
c***********************************************************************
      ! Invariants and their derivations calculation
      subroutine kInvariant(VI,xI,dIdV,ddIdV,kflagI,ndi,nshr)

      include 'aba_param.inc'

!     Data dictionary: declare variable types,definitions,& units

      ! Variables returned from the subroutine:
      real*8 xI, dIdV(6), ddIdV(6,6) ! Invariant and it's derivations

      ! Variables to be provided to the subroutine:
      real*8 V(6),VI(6)    ! Principal values

      ! Other variables:
      integer kflagI ! 1st invariant (kflagI=1), 2nd invariant (kflagI=2), 3rd invariant (kflagI=3)
      integer k, j   ! Index

      ! Initiation
      V=0.d0
      V(1:ndi)=VI(1:ndi)
      V(4:(3+nshr))=VI((ndi+1):(ndi+nshr))

      if (kflagI .eq. 1) then ! 1st Invariant

       xI=V(1)+V(2)+V(3)

!     Derivation 
       dIdV=0.d0
       dIdV(1:3)=1.d0

!     2nd derivation 
       ddIdV=0.d0

      elseif (kflagI .eq. 2) then ! 2nd Invariant

       xI=V(1)*V(2)+V(1)*V(3)+V(2)*V(3)-25.d-2*V(4)**2-25.d-2*V(5)**2-
     1 25.d-2*V(6)**2

!     Derivation 
       dIdV(1)=V(2)+V(3)
       dIdV(2)=V(1)+V(3)
       dIdV(3)=V(1)+V(2)
       dIdV(4)=-5.d-1*V(4)
       dIdV(5)=-5.d-1*V(5)
       dIdV(6)=-5.d-1*V(6)

!     2nd derivation
       ddIdV=0.d0
       do k=1,3
        do j=1,3
         if (k .ne. j) ddIdV(k,j)=1.d0
        end do
       end do

       ddIdV(4,4)=-5.d-1
       ddIdV(5,5)=-5.d-1
       ddIdV(6,6)=-5.d-1

      elseif (kflagI .eq. 3) then ! 3rd Invariant

       xI=V(1)*V(2)*V(3)+25.d-2*V(4)*V(5)*V(6)-25.d-2*V(4)**2*V(3)-
     1 25.d-2*V(6)**2*V(1)-25.d-2*V(5)**2*V(2)

!     Derivation 
       dIdV(1)=V(2)*V(3)-25.d-2*V(6)**2
       dIdV(2)=V(1)*V(3)-25.d-2*V(5)**2
       dIdV(3)=V(1)*V(2)-25.d-2*V(4)**2
       dIdV(4)=25.d-2*V(5)*V(6)-5.d-1*V(4)*V(3)
       dIdV(5)=25.d-2*V(4)*V(6)-5.d-1*V(5)*V(2)
       dIdV(6)=25.d-2*V(4)*V(5)-5.d-1*V(6)*V(1)


!     2nd derivation 
       ddIdV=0.d0
       ddIdV(1,2)=V(3)
       ddIdV(1,3)=V(2)
       ddIdV(1,6)=-5.d-1*V(6)
       ddIdV(2,3)=V(1)
       ddIdV(2,5)=-5.d-1*V(5)
       ddIdV(3,4)=-5.d-1*V(4)
       ddIdV(4,5)=25.d-2*V(6)
       ddIdV(4,6)=25.d-2*V(5)
       ddIdV(5,6)=25.d-2*V(4)

       ddIdV=ddIdV+transpose(ddIdV)

       ddIdV(4,4)=-5.d-1*V(3)
       ddIdV(5,5)=-5.d-1*V(2)
       ddIdV(6,6)=-5.d-1*V(1)

      end if

      end subroutine kInvariant

c***********************************************************************
      ! Rotating a 4th order tensor
      subroutine kRotation4thOrder(DNew,DOld,A)

      include 'aba_param.inc'

!     Data dictionary: declare variable types,definitions,& units

      ! Variables returned from the subroutine:
      real*8 DNew(6,6) ! 4th-order tensor to be rotated

      ! Variables to be provided to the subroutine:
      real*8 DOld(6,6) ! original 4th-order tensor
      real*8 A(3,3)    ! Rotation Matrix

      ! Other variables:
      integer p2u(3,3), m2i(6), m2j(6) ! Converting index from 4th-order tensor to Voigt notation and vice versa
      integer i,j,k,l,m,n,p,q,r,s,u,v  ! Index

!     Initiation
      DNew=0.d0

!     Index converting from 4th-order tensor to Voigt notation
      p2u=0
      p2u(1,2)=4
      p2u(1,3)=5
      p2u(2,3)=6

      p2u=p2u+TRANSPOSE(p2u)

      p2u(1,1)=1
      p2u(2,2)=2
      p2u(3,3)=3

!     Index converting from Voigt notation to 4th-order tensor 
      m2i=(/1, 2, 3, 1, 1, 2/)
      m2j=(/1, 2, 3, 2, 3, 3/)

!     Rotating 4th-order tensor based on rotation matrix
      do m=1,6
       do n=m,6

        i=m2i(m)
        j=m2j(m)

        k=m2i(n)
        l=m2j(n)

        DNew(m,n)=0.d0

        do p=1,3
         do q=1,3
          do r=1,3
           do s=1,3
            u=p2u(p,q)
            v=p2u(r,s)
            DNew(m,n)=DNew(m,n)+A(i,p)*A(j,q)*A(k,r)*A(l,s)*DOld(u,v)
           end do
          end do
         end do
        end do
        DNew(n,m)=DNew(m,n)
       end do
      end do

      end subroutine kRotation4thOrder

c***********************************************************************

