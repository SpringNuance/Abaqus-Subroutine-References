! UMAT and HETVAL subroutines for implementing phase field fracture
! The code is distributed under a BSD license     
      
! If using this code for research or industrial purposes, please cite:
! Navidtehrani, Y.; Betegón, C.; Martínez-Pañeda, E. A Unified Abaqus
! Implementation of the Phase Field Fracture Method Using Only a User
! Material Subroutine. Materials 2021, 14, 1913. 
! https://doi.org/10.3390/ma14081913
      
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
     1 drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,predef,dpred,
     2 cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      include 'aba_param.inc'

      character*8 cmname
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1 ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),time(2),
     2 predef(1),dpred(1),props(nprops),coords(3),drot(3,3),dfgrd0(3,3),
     3 dfgrd1(3,3),jstep(4)
     
      dimension PS(3), Edev(ntens), AN(3,3)
      
      real*8, parameter :: pi = 3.1415926535897932384626433832d0

!     Initialization
      ddsdde=0.d0
      Hmin=0.d0
      E=props(1) ! Young's modulus
      xnu=props(2) ! Poisson's ratio
      xl=props(3) ! Phase field length scale
      Gc=props(4) ! Toughness
      xk=1.d-07
      kflagS=int(props(5)) ! Solution flag (0: monolithic, 1: staggered)
      kflagC=int(props(6)) ! Model (0:AT2, 1:AT1, 2:PF-CZM/linear, 3:PF-CZM/exp])
      kflagD=int(props(7)) ! Split (0: No split, 1: Amor, 2: Miehe)
      kflagH=int(props(8)) ! Split linear momentum (1:Hybrid/no split, 2:Anisotropic [for kflagD=2])
      if (kflagC.eq.2.or.kflagC.eq.3) then
       ft=props(9) ! Tensile strength
       a=4.d0*E*Gc/(pi*xl*ft**2)
      end if 
      phi=temp+dtemp     
      psit=statev(1)
      
!     Degradation function
      call KDegFun(kflagC,phi,xk,a,g,dg,ddg)  

!     Build stiffness matrix
      call KStiffMat(E,xnu,eg,elam,ndi,nshr,ntens,bk,stran,ddsdde,
     1 kflagH,g)  
      
!     Update strains and stresses      
      stran=stran+dstran
      stress=matmul(ddsdde,stran)
      
!     Compute the driving force for fracture
      call KDriFor(E,xnu,eg,elam,bk,Gc,xl,ft,ndi,nshr,ntens,stran,
     1 stress,kflagC,kflagD,psit,g,H,Hmin)
      
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
      subroutine KStiffMat(E,xnu,eg,elam,ndi,nshr,ntens,bk,stran,ddsdde,
     1 kflagH,g)

      include 'aba_param.inc'

      dimension stress(ntens),ddsdde(ntens,ntens),stran(ntens),PS(3),
     1 AN(3,3)

      eg=E/(1.d0+xnu)/2.d0
      elam=(E/(1.d0-2.d0*xnu)-2.d0*eg)/3.d0
      bk=E/(1.d0-2.d0*xnu)/3.d0
      
      if ((kflagH.eq.1).or.(all(stran.eq.0.d0))) then
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
      elseif (kflagH .eq. 2) then
       call SPRIND(stran,PS,AN,2,ndi,nshr)
       call kddsdde(PS,AN,ddsdde,ntens,g,eg,elam)
      endif

      end
c***********************************************************************

!     Calculating 4th-order anisotropic elasticity matrix (Miehe's split)
      subroutine kddsdde(PS,AN,ddsdde,ntens,deg,eg,elam)

      include 'aba_param.inc'

      dimension kkddsdde(6,6),PS(3),ddsdde(ntens,ntens),AN(3,3)

      kkddsdde=0.d0

!     trace of strain
      TRP=PS(1)+PS(2)+PS(3)
      if (TRP.ge.0.d0) then
       kkddsdde(1:3,1:3)=elam*deg
      else
       kkddsdde(1:3,1:3)=elam
      endif

      do ki=1,3
       do kj=1,3
        do kk=1,3
         do kl=1,3
          if (ki.le.kj .and. kk.le.kl) then
           if (ki.eq.kj) then
            k1=ki
           else
            select case(10*ki+kj)
            case (12)
             k1=4
            case (23)
             k1=6
            case (13)
             k1=5
            end select
           endif
           if (kk.eq.kl) then
            k2=kk
           else
            select case(10*kk+kl)
            case (12)
             k2=4
            case (23)
             k2=6
            case (13)
             k2=5
            end select
           endif

           xKPp=0.d0
           xKPn=0.d0
           do ka=1,3
            do kb=1,3
             if (ka.eq.kb) then
              if (PS(ka) .gt. 0.d0) then
               kH=1.d0
              else
               kH=0.d0
              endif
              xkNP1=AN(ka,ki)*AN(ka,kj)*AN(kb,kk)*AN(kb,kl)
              xKPp=xKPp+kH*xkNP1
              xKPn=xKPn+(1.d0-kH)*xkNP1
             else
              xkNP1=AN(ka,ki)*AN(kb,kj)*(AN(ka,kk)*AN(kb,kl)+AN(kb,
     1        kk)*AN(ka,kl))
              if (PS(ka)==PS(kb)) then
               if(PS(ka).gt.0.d0) then
                xKPp=xKPp+5.d-1*xkNP1
               else
                xKPn=xKPn+5.d-1*xkNP1
               endif
              else
               xKPp=xKPp+5.d-1*(((PS(ka)+abs(PS(ka)))/2-(PS(kb)+
     1         abs(PS(kb)))/2)/(PS(ka)-PS(kb)))*xkNP1

               xKPn=xKPn+5.d-1*(((PS(ka)-abs(PS(ka)))/2-(PS(kb)-
     1         abs(PS(kb)))/2)/(PS(ka)-PS(kb)))*xkNP1
              endif

             endif
            end do
           end do
           kkddsdde(k1,k2)=kkddsdde(k1,k2)+2.d0*eg*(deg*xKPp+xKPn)
          endif
         end do
        end do
       end do
      end do

      ddsdde=kkddsdde(1:ntens,1:ntens)
      return
      end
c***********************************************************************

!     Compute the driving force for fracture
      subroutine KDriFor(E,xnu,eg,elam,bk,Gc,xl,ft,ndi,nshr,ntens,stran,
     1 stress,kflagC,kflagD,psit,g,H,Hmin)

      include 'aba_param.inc'

      dimension stress(ntens),ddsdde(ntens,ntens),stran(ntens),
     1 Edev(ntens),PS(3),AN(3,3)

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