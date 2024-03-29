! User material subroutine for power law conventional plasticity      
! The code is distributed under a BSD license     
      
! If using this code for research or industrial purposes, please cite:
! E. Martínez-Pañeda, S. Fuentes-Alonso, C. Betegón. 
! Gradient-enhanced statistical analysis of cleavage fracture
! European Journal of Mechanics - A/Solids 77, 103785 (2019)
! doi: 10.1016/j.euromechsol.2019.05.002 
      
! Emilio Martínez-Pañeda (mail@empaneda.com)
! Imperial College London
      
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1 drplde,drpldt,stran,dstran,time,dtime,temp2,dtemp,predef,dpred,
     2 cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      include 'aba_param.inc' 
      !implicit real(a-h o-z)

      character*8 cmname
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1 ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),
     2 time(2),predef(1),dpred(1),props(nprops),coords(3),drot(3,3),
     3 dfgrd0(3,3),dfgrd1(3,3),jstep(4)
      
      dimension eelas(ntens),eplas(ntens),flow(ntens),olds(ntens),
     + oldpl(ntens)
      
      parameter(toler=1.d-6,newton=20)

!     Initialization
      ddsdde=0.d0
      E=props(1) ! Young's modulus
      xnu=props(2) ! Poisson's ratio
      Sy=props(3) ! Yield stress
      xn=props(4) ! Strain hardening exponent
      
      call rotsig(statev(1),drot,eelas,2,ndi,nshr)
      call rotsig(statev(ntens+1),drot,eplas,2,ndi,nshr)
      eqplas=statev(1+2*ntens)
      olds=stress
      oldpl=eplas

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
       deqpl=0.d0
       Et=E*xn*(1.d0+E*eqplas/Sy)**(xn-1)
       do kewton=1,newton
        rhs=Smises-(3.d0*eg)*deqpl-Sf
        deqpl=deqpl+rhs/((3.d0*eg)+Et)
        Sf=Sy*(1.d0+E*(eqplas+deqpl)/Sy)**xn
        Et=E*xn*(1.d0+E*(eqplas+deqpl)/Sy)**(xn-1)
        if(abs(rhs).lt.toler*Sy) exit
       end do
       if (kewton.eq.newton) write(7,*)'WARNING: plasticity loop failed'

! update stresses and strains
       stress(1:3)=flow(1:3)*Sf+Sh
       eplas(1:3)=eplas(1:3)+3.d0/2.d0*flow(1:3)*deqpl
       eelas(1:3)=eelas(1:3)-3.d0/2.d0*flow(1:3)*deqpl
       stress(4:ntens)=flow(4:ntens)*Sf
       eplas(4:ntens)=eplas(4:ntens)+3.d0*flow(4:ntens)*deqpl
       eelas(4:ntens)=eelas(4:ntens)-3.d0*flow(4:ntens)*deqpl
       eqplas=eqplas+deqpl

!    Calculate the plastic strain energy density
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
      
!    Store strains in state variable array
      statev(1:ntens)=eelas
      statev((ntens+1):2*ntens)=eplas
      statev(1+2*ntens)=eqplas

      return
      end