! User element subroutine for the phase field model for fracture
! Linear elasticity coupled with piezoresistive effect, suitable for
! 3D C3D8
! Imperial College London
! Leonel Quinteros - Emilio Martinez-Pa?eda - Enrique Garc?a-Mac?as
      module kvisual
      implicit none
      real*8 UserVar(8,18,700000) !C3D8
      ! 18 variables =
      ! 6 Strains
      ! 6 Stresses
      ! 1 Phasefield
      ! 1 History field
      ! 1 Electrical potential
      ! 3 Electrical displacement
      integer nelem,kincK,kkiter
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

      parameter(ndim=3,ntens=6,ninpt=8,nsvint=18,ndof=5) !C3D8
!     ndim: Number of dimesions 3
!     ntens: Number of stresses and strains
!     ninpt: Number of integration points
!     nsinvt: Number variables
!     ndof: Number of degrees of freedom

      dimension wght(ninpt),dN(nnode,1),dNdz(ndim,nnode),dstran(ntens),
     1 dNdx(ndim,nnode),b(ntens,nnode*ndim),ddsdde(ntens,ntens),
     2 Edev(ntens),stress(ntens),stran(ntens),statevLocal(nsvint),
     3 cond(ndim,ndim),curr(ndim),dcurr(ndim),curr1(ndim)

!     Initialising: 
!     rhs: right hand side vector
!     amatrx : stifness matrix
!     wght : weights

      do k1=1,ndofel
       rhs(k1,1)=0.d0
      end do
      amatrx=0.d0
!     The weights of the integration points
	    wght=1.d0

!     Find number of elements
      if (dtime.eq.0.d0) then
       if (jelem.eq.1) then
        kincK=0
        kkiter=1
        nelem=jelem
       else
        if (jelem.gt.nelem) nelem=jelem
       endif
      endif

!     reading parameters
      xlc=props(3) ! Length scale
      Gc=props(4)  ! Critical energy release rate
      zeta=props(5) ! Viscous disipation
      kflag=props(6) ! Viscous Flag
      kflag2=props(7)  ! Monolithic or staggered
      xk=1.d-07
      eg=props(1)/(1.d0+props(2))/2.d0
      bk=props(1)/(1.d0-2.d0*props(2))/3.d0
      
!     viscous dissipation and iteration counter
      if (kflag.eq.1) then
       if (jelem.eq.1) then
        if (kinc.ne.kincK) then
         kincK=kinc
         kkiter=1
        else
         kkiter=kkiter+1
        endif
       endif
      endif

!     Iteration over the integration points of an element
      do kintk=1,ninpt
!     evaluate shape functions and derivatives
       call kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)  !Shape Funcs
       call kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd) !Derv
       call kbmatrix(dNdx,ntens,nnode,ndim,b)
       dvol=wght(kintk)*djac ! Volume differential

!     Calculate phase field and incremental strains from nodal values
       phi=0.d0
       dphi=0.d0
       volt=0.d0

!     Calculate nodal phi, volt and differential phi
       do inod=1,nnode
        phi=phi+dN(inod,1)*u(ndim*nnode+inod)
        dphi=dphi+dN(inod,1)*du(ndim*nnode+inod,1)
        volt=volt+dN(inod,1)*u(ndim*nnode+nnode+inod)
       end do

!     If phi is greater than 1, phi=1
       if (phi.gt.1.d0) phi=1.d0
!     Incremental strain
       dstran=matmul(b,du(1:ndim*nnode,1))


!     recover and assign state variables
       call kstatevar(kintk,nsvint,svars,statevLocal,1)
       stress=statevLocal(1:ntens)
       stran=statevLocal((ntens+1):(2*ntens))
       phin=statevLocal(2*ntens+1)
       Hn=statevLocal(2*ntens+2)

!     call umat to obtain stresses and constitutive matrix (monolithic)
       if (kflag2.eq.0) then
        call kumat(props,ddsdde,stress,dstran,ntens,statevLocal)
        stran=stran+dstran
        call piezor(props,stran,ntens,cond)
        phin=phi
       end if

!     compute strain energy density
       Edev=stran
       trE=stran(1)+stran(2)+stran(3)
       Edev(1:3)=stran(1:3)-trE/3.d0
       EdevS=Edev(1)**2+Edev(2)**2+Edev(3)**2
       do k1=1,ntens-3
        EdevS=EdevS+Edev(3+k1)**2/2.d0
       end do
       trEp=0.5d0*(trE+abs(trE))
       trEn=0.5d0*(trE-abs(trE))
       psip=0.5d0*bk*trEp**2+eg*EdevS
       psin=0.5d0*bk*trEn**2

!     call umat to obtain stresses and constitutive matrix (staggered)
       if (kflag2.eq.1) then
        call kumat(props,ddsdde,stress,dstran,ntens,statevLocal)
        stran=stran+dstran
        call piezor(props,stran,ntens,cond)
        if (dtime.eq.0.d0) phin=phi
       end if

!     enforcing Karush-Kuhn-Tucker conditions
       if (psip.gt.Hn) then
        H=psip
       else
        H=Hn
       endif

!     store state variables
       statevLocal(1:ntens)=stress(1:ntens)
       statevLocal((ntens+1):(2*ntens))=stran(1:ntens)
       statevLocal(2*ntens+1)=phi
       statevLocal(2*ntens+2)=H
       statevLocal(2*ntens+3)=volt
       call kstatevar(kintk,nsvint,svars,statevLocal,0)
!     form and assemble stiffness matrix and internal force vector
        g_val1=((1.d0-phin)**2+xk)
       amatrx(1:ndim*nnode,1:ndim*nnode)=
     1 amatrx(1:ndim*nnode,1:ndim*nnode)+
     2 dvol*g_val1*matmul(matmul(transpose(b),ddsdde),b)

       rhs(1:ndim*nnode,1)=rhs(1:ndim*nnode,1)-
     1 dvol*g_val1*(matmul(transpose(b),stress))

       if (kflag.eq.1) then
        if (kkiter.le.3) zeta=0.d0
       elseif (kflag.eq.2) then
        if (dphi.ge.0.d0) zeta=0.d0
       endif
!      Phase-field
       amatrx((ndim*nnode+1):(ndim*nnode+nnode),(ndim*nnode+1):
     1 (ndim*nnode+nnode))=amatrx((ndim*nnode+1):(ndim*nnode+nnode)
     2 ,(ndim*nnode+1):(ndim*nnode+nnode))+dvol*
     3 (matmul(transpose(dNdx),dNdx)*Gc*xlc+
     4 matmul(dN,transpose(dN))*(Gc/xlc+zeta/dtime+2.d0*H))

       rhs((ndim*nnode+1):(ndim*nnode+nnode),1)=
     1 rhs((ndim*nnode+1):(ndim*nnode+nnode),1)-
     2 dvol*(matmul(transpose(dNdx),matmul(dNdx
     3 ,u((ndim*nnode+1):(ndim*nnode+nnode))))
     4 *Gc*xlc+dN(:,1)*((Gc/xlc+2.d0*H)*phi+zeta*dphi/dtime-2.d0*H))

!      Degradation function that couple the piezoresistive effect

        val1= (1-exp(-(props(12))*(1-phin)**(props(13))))
        val2=  ((1-exp(-props(12))))
        g_val=(val1/val2)+xk




!      Electrical field
       amatrx((ndim*nnode+nnode+1):(ndim*nnode+2*nnode),
     1 (ndim*nnode+nnode+1):(ndim*nnode+2*nnode))=
     2 amatrx((ndim*nnode+nnode+1):(ndim*nnode+2*nnode),
     3 (ndim*nnode+nnode+1):(ndim*nnode+2*nnode))
     4 +dvol*g_val*
     5 matmul(matmul(transpose(dNdx),cond),dNdx)

       rhs((ndim*nnode+nnode+1):(ndim*nnode+2*nnode),1)=
     1 rhs((ndim*nnode+nnode+1):(ndim*nnode+2*nnode),1)
     2 -dvol*g_val*
     3 (matmul(matmul(matmul(transpose(dNdx),cond),dNdx)
     4 ,u((ndim*nnode+nnode+1):(ndim*nnode+2*nnode))))

! output
       UserVar(kintk,1:ntens,jelem)=stress(1:ntens)
       UserVar(kintk,(ntens+1):(2*ntens+3),jelem)=
     1 statevLocal((ntens+1):(2*ntens+3))
       UserVar(kintk,(2*ntens+4):(2*ntens+6),jelem)=
     1 -g_val*matmul(cond,matmul(dNdx,
     2 u((ndim*nnode+1+nnode):(ndim*nnode+2*nnode))))
       statevLocal(2*ntens+3)=volt*g_val

      end do       ! end loop on material integration points

      RETURN
      END

!***********************************************************************
      subroutine kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)
c
      include 'aba_param.inc'
c
      parameter (gaussCoord=0.577350269d0)
      dimension dN(nnode,1),dNdz(ndim,*),coord38(3,8)

      data  coord38 /-1.d0, -1.d0, -1.d0,
     2                1.d0, -1.d0, -1.d0,
     3               -1.d0,  1.d0, -1.d0,
     4                1.d0,  1.d0, -1.d0,
     5               -1.d0, -1.d0,  1.d0,
     6                1.d0, -1.d0,  1.d0,
     7               -1.d0,  1.d0,  1.d0,
     8                1.d0,  1.d0,  1.d0/

!     determine (g,h,r)
       f=coord38(1,kintk)*gaussCoord
       g=coord38(2,kintk)*gaussCoord
       h=coord38(3,kintk)*gaussCoord

!     shape functions
       dN(1,1)=0.125d0*(1.d0-f)*(1.d0-g)*(1.d0-h)
       dN(2,1)=0.125d0*(1.d0+f)*(1.d0-g)*(1.d0-h)
       dN(3,1)=0.125d0*(1.d0+f)*(1.d0+g)*(1.d0-h)
       dN(4,1)=0.125d0*(1.d0-f)*(1.d0+g)*(1.d0-h)
       dN(5,1)=0.125d0*(1.d0-f)*(1.d0-g)*(1.d0+h)
       dN(6,1)=0.125d0*(1.d0+f)*(1.d0-g)*(1.d0+h)
       dN(7,1)=0.125d0*(1.d0+f)*(1.d0+g)*(1.d0+h)
       dN(8,1)=0.125d0*(1.d0-f)*(1.d0+g)*(1.d0+h)

!     derivative d(Ni)/d(f)
       dNdz(1,1)=-0.125d0*(1.d0-g)*(1.d0-h)
       dNdz(1,2)= 0.125d0*(1.d0-g)*(1.d0-h)
       dNdz(1,3)= 0.125d0*(1.d0+g)*(1.d0-h)
       dNdz(1,4)=-0.125d0*(1.d0+g)*(1.d0-h)
       dNdz(1,5)=-0.125d0*(1.d0-g)*(1.d0+h)
       dNdz(1,6)= 0.125d0*(1.d0-g)*(1.d0+h)
       dNdz(1,7)= 0.125d0*(1.d0+g)*(1.d0+h)
       dNdz(1,8)=-0.125d0*(1.d0+g)*(1.d0+h)

!     derivative d(Ni)/d(g)
       dNdz(2,1)=-0.125d0*(1.d0-f)*(1.d0-h)
       dNdz(2,2)=-0.125d0*(1.d0+f)*(1.d0-h)
       dNdz(2,3)= 0.125d0*(1.d0+f)*(1.d0-h)
       dNdz(2,4)= 0.125d0*(1.d0-f)*(1.d0-h)
       dNdz(2,5)=-0.125d0*(1.d0-f)*(1.d0+h)
       dNdz(2,6)=-0.125d0*(1.d0+f)*(1.d0+h)
       dNdz(2,7)= 0.125d0*(1.d0+f)*(1.d0+h)
       dNdz(2,8)= 0.125d0*(1.d0-f)*(1.d0+h)

!     derivative d(Ni)/d(h)
       dNdz(3,1)=-0.125d0*(1.d0-f)*(1.d0-g)
       dNdz(3,2)=-0.125d0*(1.d0+f)*(1.d0-g)
       dNdz(3,3)=-0.125d0*(1.d0+f)*(1.d0+g)
       dNdz(3,4)=-0.125d0*(1.d0-f)*(1.d0+g)
       dNdz(3,5)= 0.125d0*(1.d0-f)*(1.d0-g)
       dNdz(3,6)= 0.125d0*(1.d0+f)*(1.d0-g)
       dNdz(3,7)= 0.125d0*(1.d0+f)*(1.d0+g)
       dNdz(3,8)= 0.125d0*(1.d0-f)*(1.d0+g)


      return
      end
!***********************************************************************
      subroutine piezor(props,stran,ntens,cond)
      include 'aba_param.inc'
      double precision rho0,pi11,pi12,pi44,pi(6,6),drho(6),rho(3,3),
     1 rhov0(6),rhof(6),detrh,cofac(3,3),stran(ntens),cond(3,3),props(*)
!      initialization of conductivity, resistivity and PI tensors
      cond = 0.d0
      pi=0.d0
      rhov0=0.d0
!      calculate the resistivity
      rho0=1/props(8)
      rhov0(1)=rho0
      rhov0(2)=rho0
      rhov0(3)=rho0
!      Extract the piezoresistive coefficients
      pi11=props(9)
      pi12=props(10)
      pi44=props(11)

!      Assign the piezoresistive coefficients to the matrix
      pi(1,1)=pi11
      pi(2,2)=pi11
      pi(3,3)=pi11
      pi(1,2)=pi12
      pi(2,1)=pi12
      pi(1,3)=pi12
      pi(3,1)=pi12
      pi(2,3)=pi12
      pi(3,2)=pi12
      pi(4,4)=pi44
      pi(5,5)=pi44
      pi(6,6)=pi44
!      Multiply the piezoresistive tensor with the strain, to then
!      multiply it with the electrical resistivity

      drho=rho0*matmul(pi,stran)
!      The isotropic coeficient is added to the piezoresistivity
!      effect
      rhof=rhov0+drho
!      The inverse of the resistivity tensor is calculated to
!      obtain the conductivity tensor
      rho(1,1)=rhof(1)
      rho(2,2)=rhof(2)
      rho(3,3)=rhof(3)
      rho(1,2)=rhof(6)
      rho(1,3)=rhof(5)
      rho(2,3)=rhof(4)
      rho(2,1)=rhof(6)
      rho(3,1)=rhof(5)
      rho(3,2)=rhof(4)
      detrh =  rho(1,1)*rho(2,2)*rho(3,3)
     1      - rho(1,1)*rho(2,3)*rho(3,2)
     2      - rho(1,2)*rho(2,1)*rho(3,3)
     3      + rho(1,2)*rho(2,3)*rho(3,1)
     4      + rho(1,3)*rho(2,1)*rho(3,2)
     5      - rho(1,3)*rho(2,2)*rho(3,1)


      cofac(1,1) = +(rho(2,2)*rho(3,3)-rho(2,3)*rho(3,2))
      cofac(1,2) = -(rho(2,1)*rho(3,3)-rho(2,3)*rho(3,1))
      cofac(1,3) = +(rho(2,1)*rho(3,2)-rho(2,2)*rho(3,1))
      cofac(2,1) = -(rho(1,2)*rho(3,3)-rho(1,3)*rho(3,2))
      cofac(2,2) = +(rho(1,1)*rho(3,3)-rho(1,3)*rho(3,1))
      cofac(2,3) = -(rho(1,1)*rho(3,2)-rho(1,2)*rho(3,1))
      cofac(3,1) = +(rho(1,2)*rho(2,3)-rho(1,3)*rho(2,2))
      cofac(3,2) = -(rho(1,1)*rho(2,3)-rho(1,3)*rho(2,1))
      cofac(3,3) = +(rho(1,1)*rho(2,2)-rho(1,2)*rho(2,1))
      cond = transpose(cofac)/detrh
      return
      end

!***********************************************************************
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

      if (ndim.eq.3) then

       djac=xjac(1,1)*xjac(2,2)*xjac(3,3)+xjac(2,1)*xjac(3,2)*xjac(1,3)
     & +xjac(3,1)*xjac(2,3)*xjac(1,2)-xjac(3,1)*xjac(2,2)*xjac(1,3)
     & -xjac(2,1)*xjac(1,2)*xjac(3,3)-xjac(1,1)*xjac(2,3)*xjac(3,2)
       if (djac.gt.0.d0) then ! jacobian is positive - o.k.
        xjaci(1,1)=(xjac(2,2)*xjac(3,3)-xjac(2,3)*xjac(3,2))/djac
        xjaci(1,2)=(xjac(1,3)*xjac(3,2)-xjac(1,2)*xjac(3,3))/djac
        xjaci(1,3)=(xjac(1,2)*xjac(2,3)-xjac(1,3)*xjac(2,2))/djac
        xjaci(2,1)=(xjac(2,3)*xjac(3,1)-xjac(2,1)*xjac(3,3))/djac
        xjaci(2,2)=(xjac(1,1)*xjac(3,3)-xjac(1,3)*xjac(3,1))/djac
        xjaci(2,3)=(xjac(1,3)*xjac(2,1)-xjac(1,1)*xjac(2,3))/djac
        xjaci(3,1)=(xjac(2,1)*xjac(3,2)-xjac(2,2)*xjac(3,1))/djac
        xjaci(3,2)=(xjac(1,2)*xjac(3,1)-xjac(1,1)*xjac(3,2))/djac
        xjaci(3,3)=(xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1))/djac
       else ! negative or zero jacobian
        write(7,*)'WARNING: element',jelem,'has neg. Jacobian'
       endif

      else if (ndim.eq.2) then

       djac=xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1)
       if (djac.gt.0.d0) then ! jacobian is positive - o.k.
        xjaci(1,1)=xjac(2,2)/djac
        xjaci(2,2)=xjac(1,1)/djac
        xjaci(1,2)=-xjac(1,2)/djac
        xjaci(2,1)=-xjac(2,1)/djac
       else ! negative or zero jacobian
        write(7,*)'WARNING: element',jelem,'has neg. Jacobian'
       endif

      endif

      dNdx=matmul(xjaci,dNdz)

      return
      end

!***********************************************************************
      subroutine kbmatrix(dNdx,ntens,nnode,ndim,b)
!     Notation, strain tensor: e11, e22, e33, e12, e13, e23
      include 'aba_param.inc'

      dimension dNdx(ndim,nnode),b(ntens,nnode*ndim)

      b=0.d0
      do inod=1,nnode
       b(1,ndim*inod-ndim+1)=dNdx(1,inod)
       b(2,ndim*inod-ndim+2)=dNdx(2,inod)
       b(4,ndim*inod-ndim+1)=dNdx(2,inod)
       b(4,ndim*inod-ndim+2)=dNdx(1,inod)
       if (ndim.eq.3) then
        b(3,ndim*inod)=dNdx(3,inod)
        b(5,ndim*inod-2)=dNdx(3,inod)
        b(5,ndim*inod)=dNdx(1,inod)
        b(6,ndim*inod-1)=dNdx(3,inod)
        b(6,ndim*inod)=dNdx(2,inod)
       endif
      end do

      return
      end

!***********************************************************************
      subroutine kstatevar(npt,nsvint,statev,statev_ip,icopy)
!
!     Transfer data to/from element-level state variable array from/to
!     material-point level state variable array.
!
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

!***********************************************************************
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
      do k1=4,ntens
       ddsdde(k1,k1)=eg2/2.d0
      end do

      stress=stress+matmul(ddsdde,dstran)

      return
      end

!***********************************************************************
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
      statev(1:nstatv)=UserVar(npt,1:nstatv,noffset)

      return
      end
