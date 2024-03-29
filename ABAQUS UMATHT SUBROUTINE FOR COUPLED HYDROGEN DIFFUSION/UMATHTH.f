! User material and thermal subroutine for hydrogen transport
! The code is distributed under a BSD license
      
! If using this code for research or industrial purposes, please cite:
! R. Fernandez-Sousa, C. Betegon, E. Martinez-Paneda. Analysis of the 
! influence of microstructural traps on hydrogen assisted fatigue.
! Acta Materialia 199, pp. 253-263 (2020). 
! doi: 10.1016/j.actamat.2020.08.030
      
! Emilio Martinez-Paneda (e.martinez-paneda@imperial.ac.uk)
! Imperial College London
      
      module ktransfer
      implicit none
      real*8 coorT(50000,4,2),grad(50000,8),ShT(50000,4)
      save
      end module   

!***********************************************************************
      subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)
      use ktransfer
      include 'aba_param.inc' !implicit real(a-h o-z)
      dimension time(2)
      
      if (lop.eq.0) then ! start of the analysis
       coorT=0.d0
       grad=0.d0
       ShT=0.d0
      end if
      
      return
      end
      
!***********************************************************************      
      subroutine umatht(u,dudt,dudg,flux,dfdt,dfdg,statev,temp,dtemp,
     1 dtemdx,time,dtime,predef,dpred,cmname,ntgrd,nstatv,props,nprops,
     2 coords,pnewdt,noel,npt,layer,kspt,kstep,kinc)

      use ktransfer
      include 'aba_param.inc'

      character*80 cmname
      dimension dudg(ntgrd),flux(ntgrd),dfdt(ntgrd),dfdg(ntgrd,ntgrd),
     1 statev(nstatv),dtemdx(ntgrd),time(2),predef(1),dpred(1),
     2 props(nprops),coords(3)

      dimension Wb(nprops/2-1),xK(nprops/2-1),xNt(nprops/2-1),sig(ntgrd)

      ! Step-1: Read input data & Initialize
      dfdg=0.d0
      ntens=ntgrd*2
      cL=temp+dtemp

      D=props(1)
      kflag=props(2)
      xNl=5.1d20 ! [sites/mm^3]
      Vh=2000.d0 ! [mm^3/mol]
      R=8314.5d0 ! [N*mm/(mol*K)]
      T=300.d0   ! [K]
      b=0.2725d-6 ! bcc [mm]
      sig(1)=grad(noel,2*npt-1)
      sig(2)=grad(noel,2*npt)
      ntraps=nprops/2-1
      
      do k1=1,ntraps
       Wb(k1)=props(2+2*k1-1)
       xNt(k1)=props(2+2*k1)
       xK(k1)=exp(Wb(k1)/(R*T))
      end do
      
      du2=0.d0
      if (kflag.eq.0) then ! Lattice H only
       xNt(1)=0.d0
      elseif (kflag.eq.1) then ! No dislocations: Wb(1) and xNt(1) irrelevant
       xNt(1)=0.d0
      elseif (kflag.eq.2) then ! Kumnick & Johnson / Sofronis & McMeeking (in sites/mm^3)
       xNt(1)=10.d0**(23.26d0-2.33d0*exp(-5.5d0*statev(1+2*ntens)))/1e9
      elseif (kflag.eq.3) then ! Kumnick & Johnson / Krom et al. 
       xNt(1)=10.d0**(23.26d0-2.33d0*exp(-5.5d0*statev(1+2*ntens)))/1e9
       du2=(xK(1)*cL/(xK(1)*cL+xNl))*29.5d0
     & *dexp(-5.5d0*statev(1+2*ntens))*xNt(1)*statev(2+2*ntens)
      elseif (kflag.eq.4) then ! Gilman / Dadfarnia et al.
        if (statev(1+2*ntens).lt.0.5) then
         xNt(1)=(1.d10+statev(1+2*ntens)*2.d16)/(b*1e6)
         du2=(xK(1)*cL/(xK(1)*cL+xNl))*(statev(2+2*ntens)*2e16)/(b*1e6)
        elseif (statev(1+2*ntens).ge.0.5) then
         xNt(1)=(1e16)/(b*1e6)
        endif
      elseif (kflag.eq.5) then ! Taylor / Fernandez-Sousa et al.
       xNt(1)=statev(3+2*ntens)/b
      endif
      
      dudt2=0.d0
      do k1=1,ntraps
       dudt2=dudt2+xNt(k1)*xK(k1)*xNl/((xK(k1)*cL+xNl)**2.d0)
      end do 
      dudt=1.d0+dudt2	   
      u=u+dudt*dtemp+du2
      do i=1,ntgrd
       dudg(i)=0.d0
       flux(i)=-D*dtemdx(i)+D*cL*Vh*sig(i)/(R*T)
       dfdt(i)=D*Vh*sig(i)/(R*T)
       dfdg(i,i)=-D
      end do
      
!     store the concentration in each trap, in all traps and in traps and lattice
      id=3+2*ntens
      statev(ntraps+1+id)=0
      do k1=1,ntraps      
       statev(k1+id)=xNt(k1)*xK(k1)*cL/(xK(k1)*cL+xNl)
       statev(ntraps+1+id)=statev(ntraps+1+id)+statev(k1+id)
      end do
      statev(ntraps+2+id)=cL+statev(ntraps+1+id)
	  
      return
      end
      
! User material subroutine for power law conventional plasticity      
! The code is distributed under a BSD license     
      
! If using this code for research or industrial purposes, please cite:
! E. Martinez-Paneda, S. Fuentes-Alonso, C. Betegon. 
! Gradient-enhanced statistical analysis of cleavage fracture
! European Journal of Mechanics - A/Solids 77, 103785 (2019)
! doi: 10.1016/j.euromechsol.2019.05.002 
      
! Emilio Martinez-Paneda (mail@empaneda.com)
! Imperial College London
      
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1 drplde,drpldt,stran,dstran,time,dtime,temp2,dtemp,predef,dpred,
     2 cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      use ktransfer
      include 'aba_param.inc'

      character*8 cmname
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1 ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),
     2 time(2),predef(1),dpred(1),props(nprops),coords(3),drot(3,3),
     3 dfgrd0(3,3),dfgrd1(3,3),jstep(4)
      
      dimension eelas(ntens),eplas(ntens),flow(ntens),olds(ntens),
     + oldpl(ntens),deriv(2,4),xjacm(2,2),xjaci(2,2)
           
      parameter(toler=1.d-6,newton=20)

!     Initialization
      ddsdde=0.d0
      deqpl=0.d0
      E=props(1) ! Young's modulus
      xnu=props(2) ! Poisson's ratio
      Sy=props(3) ! Yield stress
      xn=props(4) ! Strain hardening exponent
      
      call rotsig(statev(1),drot,eelas,2,ndi,nshr)
      call rotsig(statev(ntens+1),drot,eplas,2,ndi,nshr)
      eqplas=statev(1+2*ntens)
      olds=stress
      oldpl=eplas

!     Compute the gradient of the hydrostatic stress    
      if (npt==1 .and. time(1).gt.0) then
       do k1=1,4 ! Hard coded for 4 integration points
        if (k1==1) then
         s=-1.d0
         t=-1.d0
        elseif (k1==2) then 
         s=1.d0
         t=-1.d0
        elseif (k1==3) then
         s=-1.d0
         t=1.d0  
        elseif (k1==4) then      
         s=1.d0
         t=1.d0   
        end if

        deriv(1,1)=-(1.d0/4.0)*(1-t)
        deriv(1,2)=(1.d0/4.0)*(1-t)
        deriv(1,3)=-(1.d0/4.0)*(1+t)
        deriv(1,4)=(1.d0/4.0)*(1+t)
        deriv(2,1)=-(1.d0/4.0)*(1-s)
        deriv(2,2)=-(1.d0/4.0)*(1+s)
        deriv(2,3)=(1.d0/4.0)*(1-s)
        deriv(2,4)=(1.d0/4.0)*(1+s)

        xjacm(1,1)=deriv(1,1)*coorT(noel,1,1)+deriv(1,2)*coorT(noel,2,1)
     1 +deriv(1,3)*coorT(noel,3,1)+deriv(1,4)*coorT(noel,4,1)
    
        xjacm(1,2)=deriv(1,1)*coorT(noel,1,2)+deriv(1,2)*coorT(noel,2,2)
     1 +deriv(1,3)*coorT(noel,3,2)+deriv(1,4)*coorT(noel,4,2)
     
        xjacm(2,1)=deriv(2,1)*coorT(noel,1,1)+deriv(2,2)*coorT(noel,2,1)
     1 +deriv(2,3)*coorT(noel,3,1)+deriv(2,4)*coorT(noel,4,1)
      
        xjacm(2,2)=deriv(2,1)*coorT(noel,1,2)+deriv(2,2)*coorT(noel,2,2)
     1 +deriv(2,3)*coorT(noel,3,2)+deriv(2,4)*coorT(noel,4,2)

        djacb=xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1) 
      
        xjaci(1,1)=xjacm(2,2)/djacb 
        xjaci(1,2)=-xjacm(1,2)/djacb  
        xjaci(2,1)=-xjacm(2,1)/djacb   
        xjaci(2,2)=xjacm(1,1)/djacb

        a1=xjaci(1,1)*deriv(1,1)+xjaci(1,2)*deriv(2,1) 
        a2=xjaci(1,1)*deriv(1,2)+xjaci(1,2)*deriv(2,2) 
        a3=xjaci(1,1)*deriv(1,3)+xjaci(1,2)*deriv(2,3) 
        a4=xjaci(1,1)*deriv(1,4)+xjaci(1,2)*deriv(2,4) 
        b1=xjaci(2,1)*deriv(1,1)+xjaci(2,2)*deriv(2,1) 
        b2=xjaci(2,1)*deriv(1,2)+xjaci(2,2)*deriv(2,2) 
        b3=xjaci(2,1)*deriv(1,3)+xjaci(2,2)*deriv(2,3) 
        b4=xjaci(2,1)*deriv(1,4)+xjaci(2,2)*deriv(2,4)  
      
        grad(noel,2*k1-1)=a1*ShT(noel,1)+a2*ShT(noel,2)+a3*ShT(noel,3)
     1 +a4*ShT(noel,4)
        grad(noel,2*k1)=b1*ShT(noel,1)+b2*ShT(noel,2)+b3*ShT(noel,3)
     1 +b4*ShT(noel,4)
       end do 
      end if     
      
!     Build elastic stiffness matrix
      eg=E/(1.d0+xnu)/2.d0
      elam=(E/(1.d0-2.d0*xnu)-2.d0*eg)/3.d0
      
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
      stress=stress+matmul(ddsdde,dstran)
      eelas=eelas+dstran
      
!     Calculate equivalent von Mises stress
      Smises=(stress(1)-stress(2))**2+(stress(2)-stress(3))**2
     1 +(stress(3)-stress(1))**2
      do i=4,ntens
       Smises=Smises+6.d0*stress(i)**2
      end do
      Smises=sqrt(Smises/2.d0)	 
      
!     Get yield stress from the specified hardening curve
      Sf=Sy*(1.d0+E*eqplas/Sy)**xn
      
!     Determine if active yielding
      if (Smises.gt.(1.d0+toler)*Sf) then

!     Calculate the flow direction
       Sh=(stress(1)+stress(2)+stress(3))/3.d0
       flow(1:3)=(stress(1:3)-Sh)/Smises
       flow(4:ntens)=stress(4:ntens)/Smises
       
!     Solve for Smises and deqpl using Newton's method
       Et=E*xn*(1.d0+E*eqplas/Sy)**(xn-1)
       do kewton=1,newton
        rhs=Smises-(3.d0*eg)*deqpl-Sf
        deqpl=deqpl+rhs/((3.d0*eg)+Et)
        Sf=Sy*(1.d0+E*(eqplas+deqpl)/Sy)**xn
        Et=E*xn*(1.d0+E*(eqplas+deqpl)/Sy)**(xn-1)
        if(abs(rhs).lt.toler*Sy) exit
       end do
       if (kewton.eq.newton) write(7,*)'WARNING: plasticity loop failed'

!     Update stresses and strains
       stress(1:3)=flow(1:3)*Sf+Sh
       eplas(1:3)=eplas(1:3)+3.d0/2.d0*flow(1:3)*deqpl
       eelas(1:3)=eelas(1:3)-3.d0/2.d0*flow(1:3)*deqpl
       stress(4:ntens)=flow(4:ntens)*Sf
       eplas(4:ntens)=eplas(4:ntens)+3.d0*flow(4:ntens)*deqpl
       eelas(4:ntens)=eelas(4:ntens)-3.d0*flow(4:ntens)*deqpl
       eqplas=eqplas+deqpl

!     Calculate the plastic strain energy density
       do i=1,ntens
        spd=spd+(stress(i)+olds(i))*(eplas(i)-oldpl(i))/2.d0
       end do
      
!     Formulate the jacobian (material tangent)   
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
      
!     Store strains in state variable array
      statev(1:ntens)=eelas
      statev((ntens+1):2*ntens)=eplas
      statev(1+2*ntens)=eqplas
      statev(2+2*ntens)=deqpl
      
!     Statistically stored dislocations (bcc)
      xm=2.9d0
      b=0.2725d-6
      ssd=((Sy*(E/Sy)**xn*(eqplas+Sy/E)**xn)/(xm*0.5d0*eg*b))**2 
      statev(3+2*ntens)=ssd
      
!     Store information for computing the Sh gradient     
      ShT(noel,npt)=(stress(1)+stress(2)+stress(3))/3.d0      
      coorT(noel,npt,1)=coords(1)
      coorT(noel,npt,2)=coords(2)

      return
      end
      
! User subroutine for prescribing a remote amplitude load
! The code is distributed under a BSD license     
      
! If using this code for research or industrial purposes, please cite:
! S. del Busto, C. Betegón, E. Martínez-Pañeda. A cohesive zone 
! framework for environmentally assisted fatigue. Engineering Fracture
! Mechanics 185: 210-226 (2017). doi:10.1016/j.engfracmech.2017.05.021     
      
! Emilio Martínez-Pañeda (mail@empaneda.com)
! Technical University of Denmark
      
      subroutine  disp(u,kstep,kinc,time,node,noel,jdof,coords)

      include 'aba_param.inc'

      dimension u(3),time(2),coords(3)
      
! Introduce manually        
      xk=30*dsqrt(1000.d0)*time(2) ! MPa*m^1/2 to MPa*mm^1/2
      e=201880.d0
      xnu=0.3d0
      pi=3.14159265359d0    

!user coding to define U
      r=dsqrt(coords(1)**2+coords(2)**2)
      theta=atan2(coords(2),coords(1))
      if (jdof.eq.1) then
       u(1)=(1+xnu)*(xk/e)*sqrt(r/(2*pi))*(3-4*xnu-cos(theta))*
     & cos(theta/2)
      elseif (jdof.eq.2) then
       u(1)=(1+xnu)*(xk/e)*sqrt(r/(2*pi))*(3-4*xnu-cos(theta))*
     & sin(theta/2)
      endif  

      RETURN
      END       