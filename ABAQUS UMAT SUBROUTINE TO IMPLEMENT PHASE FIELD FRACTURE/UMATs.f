! UMAT subroutine for implementing the phase field model for fracture
! The code is distributed under a BSD license     
      
! If using this code for research or industrial purposes, please cite:
! Navidtehrani, Y.; Betegón, C.; Martínez-Pañeda, E. A Unified Abaqus
! Implementation of the Phase Field Fracture Method Using Only a User
! Material Subroutine. Materials 2021, 14, 1913. 
! https://doi.org/10.3390/ma14081913
      
! Emilio Martinez-Paneda (e.martinez-paneda@imperial.ac.uk)
! Imperial College London
c*****************************************************************
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

!     Initialization
      ddsdde=0.d0
      E=props(1) ! Young's modulus
      xnu=props(2) ! Poisson's ratio
      xl=props(3) ! Phase field length scale
      Gc=props(4) ! Toughness
      kflag=int(props(5)) ! Solution flag (0: monolithic, 1: staggered)
      phi=temp+dtemp     
      psit=statev(1)
      g=(1.d0-phi)**2+1.d-07 ! Degradation function      

!     Build stiffness matrix
      eg=E/(1.d0+xnu)/2.d0
      elam=(E/(1.d0-2.d0*xnu)-2.d0*eg)/3.d0
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
      stran=stran+dstran
      stress=matmul(ddsdde,stran)
      
!     Compute the strain energy density
      psi=0.d0
      do i=1,ntens
       psi=psi+stress(i)*stran(i)*0.5d0
      end do
      H=max(psit,psi)
      stress=stress*g
      ddsdde=ddsdde*g 
      statev(1)=H

!     phase field      
      if (kflag.eq.1) H=psit    ! Staggered scheme
      rpl=-(phi/xl**2-2.d0*(1.d0-phi)*H/(Gc*xl))
      drpldt=-(1.d0/xl**2+2.d0*H/(Gc*xl))
      
      return
      end