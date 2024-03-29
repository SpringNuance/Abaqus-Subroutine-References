! User element subroutine for the phase field model for fracture
! including coupling with hydrogen diffusion (mass transport)     
! Quadratic quadrilateral elements, explicit with eq. correction
! The code is distributed under a BSD license     
      
! If using this code for research or industrial purposes, please cite:
! E. Martínez-Pañeda, A. Golahmar and C.F. Niordson. 
! A phase field formulation for hydrogen assisted cracking
! Computer Methods in Applied Mechanics and Engineering 342: 742-761 
! (2018) doi: 10.1016/j.cma.2018.07.021
      
! Emilio Martínez-Pañeda (mail@empaneda.com)
! University of Cambridge

      module kvisual
      implicit none
      real*8 UserVar(70000,11,4)
      integer nelem
      save
      end module
      
      subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     1 props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,
     2 kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf,
     3 lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njpro,period)

      use kvisual
      include 'aba_param.inc' !implicit real(a-h o-z)
      
      dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),svars(*),
     1 energy(*),coords(mcrd,nnode),u(ndofel),du(mlvarx,*),v(ndofel),
     2 a(ndofel),time(2),params(*),jdltyp(mdload,*),adlmag(mdload,*),
     3 ddlmag(mdload,*),predef(2,npredf,nnode),lflags(*),jprops(*)

      parameter(ndim=2,ntens=4,ninpt=4,nsvint=10)
      
      dimension wght(ninpt),dN(1,nnode),dNdz(ndim,nnode),dNS(ninpt),
     2 dNdx(ndim,nnode),b(ntens,nnode*ndim),ddsdde(ntens,ntens),
     3 stress(ntens),stran(ntens),bC(ndim,nnode),xm(nnode,nnode),
     4 xk(nnode,nnode),BB(nnode,nnode),SHa(nnode,1),coord28(ndim,nnode),
     5 dstran(ntens),statevLocal(nsvint)
      
      data wght /1.d0, 1.d0, 1.d0, 1.d0/
        
!     initialising
      do k1=1,ndofel
       rhs(k1,1)=0.d0
      end do
      amatrx=0.d0
      
!     find number of elements          
      if (dtime.eq.0.d0) then
       if (jelem.eq.1) then
        nelem=jelem
       else
        if (jelem.gt.nelem) nelem=jelem 
       endif 
      endif      
      
!     reading parameters
      xlc=props(3)
      Gc0=props(4)
      D=props(5) ! Diffusion coefficient
      xkap=1.d-7 ! well-conditioning parameter
      Vh=2000.d0 ! Molar volume of H
      T=300.d0 ! Temperature
      R=8314.5d0 ! Gas constant

!     compute the hydrostatic stress
      SHa=0.d0
      coord28=0.d0
      coord28(1,1)=-1.d0
      coord28(2,1)=-1.d0
      coord28(1,2)=1.d0
      coord28(2,2)=-1.d0
      coord28(1,3)=1.d0
      coord28(2,3)=1.d0  
      coord28(1,4)=-1.d0
      coord28(2,4)=1.d0
      coord28(2,5)=-1.d0
      coord28(1,6)=1.d0
      coord28(2,7)=1.d0  
      coord28(1,8)=-1.d0        
      
      do inod=1,nnode
       g=dsqrt(3.d0)*coord28(1,inod)
       h=dsqrt(3.d0)*coord28(2,inod)
       dNS(1)=(1.d0-g)*(1.d0-h)/4.d0
       dNS(2)=(1.d0+g)*(1.d0-h)/4.d0
       dNS(3)=(1.d0-g)*(1.d0+h)/4.d0
       dNS(4)=(1.d0+g)*(1.d0+h)/4.d0
       do i=1,ninpt
        isvinc=(i-1)*nsvint
        SHa(inod,1)=SHa(inod,1)+dNS(i)*svars(isvinc+10)
       end do
      end do      
      
      do kintk=1,ninpt
!     evaluate shape functions and derivatives
       call kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)      
       call kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
       dvol=wght(kintk)*djac
       
!     form B-matrix
       b=0.d0
       do inod=1,nnode
        bC(1,inod)=dNdx(1,inod)
        bC(2,inod)=dNdx(2,inod)
        b(1,2*inod-1)=dNdx(1,inod)
        b(2,2*inod)=dNdx(2,inod)
        b(4,2*inod-1)=dNdx(2,inod)
        b(4,2*inod)=dNdx(1,inod)
       end do                     
       
!     compute from nodal values
       phi=0.d0
       cL=0.d0
       do inod=1,nnode
        phi=phi+dN(1,inod)*u(ndim*nnode+inod)
        cL=cL+dN(1,inod)*u((ndim+1)*nnode+inod)
       end do   
       if (phi.gt.1.d0) phi=1.d0
       
!     hydrogen contribution  
       Theta=cL*5.5d-05/(cL*5.5d-05+dexp(-3.d7/(R*T)))
       Gc=Gc0*(1.d0-0.89d0*Theta)
           
!     compute the increment of strain and recover history variables
       dstran=matmul(b,du(1:ndim*nnode,1))
       call kstatevar(kintk,nsvint,svars,statevLocal,1)
       stress=statevLocal(1:ntens)
       stran(1:ntens)=statevLocal((ntens+1):(2*ntens))
       Hn=statevLocal(2*ntens+1)
       
!     compute strain energy density from the previous increment       
       Psi=0.d0
       do k1=1,ntens
        Psi=Psi+stress(k1)*stran(k1)*0.5d0
       end do
       
!     call umat to obtain stresses and constitutive matrix 
       call kumat(props,ddsdde,stress,dstran,ntens,statevLocal)
       stran=stran+dstran
       
!     enforcing Karush-Kuhn-Tucker conditions
       if (Psi.gt.Hn) then
        H=Psi
       else
        H=Hn
       endif
       
       statevLocal(1:ntens)=stress(1:ntens)
       statevLocal((ntens+1):(2*ntens))=stran(1:ntens)
       statevLocal(2*ntens+1)=H
       statevLocal(2*ntens+2)=(stress(1)+stress(2)+stress(3))/3.d0 !SH
       
       call kstatevar(kintk,nsvint,svars,statevLocal,0)
       
       amatrx(1:16,1:16)=amatrx(1:16,1:16)+dvol*(((1.d0-phi)**2+xkap)
     1 *matmul(matmul(transpose(b),ddsdde),b))
        
       rhs(1:16,1)=rhs(1:16,1)-
     1 dvol*(matmul(transpose(b),stress)*((1.d0-phi)**2+xkap))       
           
       amatrx(17:24,17:24)=amatrx(17:24,17:24)
     1 +dvol*(matmul(transpose(dNdx),dNdx)*Gc*xlc
     2 +matmul(transpose(dN),dN)*(Gc/xlc+2.d0*H))   

       rhs(17:24,1)=rhs(17:24,1)
     1 -dvol*(matmul(transpose(dNdx),matmul(dNdx,u(17:24)))
     2 *Gc*xlc+dN(1,:)*((Gc/xlc+2.d0*H)*phi-2.d0*H))
        
       xm=matmul(transpose(dN),dN)/D
       BB=matmul(transpose(bC),bC)
       xk=BB-Vh/(R*T)*matmul(BB,matmul((SHa*((1.d0-phi)**2+xkap)),dN))
        
       amatrx(25:32,25:32)=amatrx(25:32,25:32)+dvol*(xm/dtime+xk)
        
       rhs(25:32,1)=rhs(25:32,1)-dvol*(matmul(xk,u(25:32))+
     1 matmul(xm,du(25:32,1))/dtime)
                   
! output
       UserVar(jelem,1:4,kintk)=statevLocal(1:4)*((1.d0-phi)**2+xkap)
       UserVar(jelem,5:8,kintk)=statevLocal((ntens+1):(2*ntens))
       UserVar(jelem,(2*ntens+1),kintk)=phi
       UserVar(jelem,(2*ntens+2),kintk)=statevLocal(2*ntens+2)
       UserVar(jelem,(2*ntens+3),kintk)=cL
      
      end do       ! end loop on material integration points
      
      RETURN
      END
      
      subroutine kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)
c
      include 'aba_param.inc'
c
      parameter (gaussCoord=0.577350269d0)
      dimension dN(1,nnode),dNdz(ndim,*),coord24(2,4)
      
      data  coord24 /-1.d0, -1.d0,
     2                1.d0, -1.d0,
     3               -1.d0,  1.d0,
     4                1.d0,  1.d0/      
     
!     2D 4-nodes

!     determine (g,h,r)
      g=coord24(1,kintk)*gaussCoord
      h=coord24(2,kintk)*gaussCoord

!     shape functions 
      dN(1,1)=-0.25d0*(1.d0-g)*(1.d0-h)*(1.d0+g+h)
      dN(1,2)=0.25d0*(1.d0+g)*(1.d0-h)*(g-h-1.d0)
      dN(1,3)=0.25d0*(1.d0+g)*(1.d0+h)*(g+h-1.d0)
      dN(1,4)=0.25d0*(1.d0-g)*(1.d0+h)*(h-g-1.d0)
      dN(1,5)=0.5d0*(1.d0-g*g)*(1.d0-h)
      dN(1,6)=0.5d0*(1.d0+g)*(1.d0-h*h)
      dN(1,7)=0.5d0*(1.d0-g*g)*(1.d0+h)
      dN(1,8)=0.5d0*(1.d0-g)*(1.d0-h*h)        

!     derivative d(Ni)/d(g)
      dNdz(1,1)=0.25d0*(1.d0-h)*(2.d0*g+h)
      dNdz(1,2)=0.25d0*(1.d0-h)*(2.d0*g-h)
      dNdz(1,3)=0.25d0*(1.d0+h)*(2.d0*g+h)
      dNdz(1,4)=0.25d0*(1.d0+h)*(2.d0*g-h)
      dNdz(1,5)=-g*(1.d0-h)
      dNdz(1,6)=0.5d0*(1.d0-h*h)
      dNdz(1,7)=-g*(1.d0+h)
      dNdz(1,8)=-0.5d0*(1.d0-h*h)      

!     derivative d(Ni)/d(h)
      dNdz(2,1)=0.25d0*(1.d0-g)*(g+2.d0*h)
      dNdz(2,2)=0.25d0*(1.d0+g)*(2.d0*h-g)
      dNdz(2,3)=0.25d0*(1.d0+g)*(2.d0*h+g)
      dNdz(2,4)=0.25d0*(1.d0-g)*(2.d0*h-g)
      dNdz(2,5)=-0.5d0*(1.d0-g*g) 
      dNdz(2,6)=-(1.d0+g)*h 
      dNdz(2,7)=0.5d0*(1.d0-g*g)
      dNdz(2,8)=-(1.d0-g)*h  
      
      return
      end 

      subroutine kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
!     Notation: djac - Jac determinant; xjaci - inverse of Jac matrix 
!     dNdx - shape functions derivatives w.r.t. global coordinates
      include 'aba_param.inc'

      dimension xjac(ndim,ndim),xjaci(ndim,ndim),coords(mcrd,nnode),
     1 dNdz(ndim,nnode),dNdx(ndim,nnode)

      xjac=0.d0

      do inod=1,nnode
       do idim=1,ndim
        do jdim=1,ndim
         xjac(jdim,idim)=xjac(jdim,idim)+
     1        dNdz(jdim,inod)*coords(idim,inod)      
        end do
       end do 
      end do

      djac=xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1)
      if (djac.gt.0.d0) then ! jacobian is positive - o.k.
       xjaci(1,1)=xjac(2,2)/djac
       xjaci(2,2)=xjac(1,1)/djac
       xjaci(1,2)=-xjac(1,2)/djac
       xjaci(2,1)=-xjac(2,1)/djac
      else ! negative or zero jacobian
       write(7,*)'WARNING: element',jelem,'has neg. Jacobian'
      endif
	  
      dNdx=matmul(xjaci,dNdz) 
			
      return
      end

c*****************************************************************
      subroutine kstatevar(npt,nsvint,statev,statev_ip,icopy)
c
c     Transfer data to/from element-level state variable array from/to
c     material-point level state variable array.
c
      include 'aba_param.inc'

      dimension statev(*),statev_ip(*)

      isvinc=(npt-1)*nsvint     ! integration point increment

      if (icopy.eq.1) then ! Prepare arrays for entry into umat
       do i=1,nsvint
        statev_ip(i)=statev(i+isvinc)
       enddo
      else ! Update element state variables upon return from umat
       do i=1,nsvint
        statev(i+isvinc)=statev_ip(i)
       enddo
      end if

      return
      end

c*****************************************************************
      subroutine kumat(props,ddsdde,stress,dstran,ntens,statev)
c
c     Subroutine with the material model
c
      include 'aba_param.inc' !implicit real(a-h o-z)
      
      dimension props(*),ddsdde(ntens,ntens),stress(ntens),statev(*),
     + dstran(ntens)

!     Initialization
      ddsdde=0.d0
      E=props(1) ! Young's modulus
      xnu=props(2) ! Poisson's ratio
      
!     Build stiffness matrix
      eg2=E/(1.d0+xnu)
      elam=(E/(1.d0-2.d0*xnu)-eg2)/3.d0
      
!     Update stresses
      do k1=1,3
       do k2=1,3
        ddsdde(k2,k1)=elam
       end do
       ddsdde(k1,k1)=eg2+elam
      end do
      ddsdde(4,4)=eg2/2.d0
      
      stress=stress+matmul(ddsdde,dstran)   

      return
      end

c*****************************************************************
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1 drplde,drpldt,stran,dstran,time,dtime,temp2,dtemp,predef,dpred,
     2 cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      use kvisual
      include 'aba_param.inc' !implicit real(a-h o-z)

      character*8 cmname
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1 ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),
     2 time(2),predef(1),dpred(1),props(nprops),coords(3),drot(3,3),
     3 dfgrd0(3,3),dfgrd1(3,3),jstep(4)

      ddsdde=0.0d0
      noffset=noel-nelem
      statev(1:nstatv)=UserVar(noffset,1:nstatv,npt)
     
      return
      end
