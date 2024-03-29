! UELMAT subroutine for the phase field fracture method
! The code is distributed under a BSD license     
      
! If using this code for research or industrial purposes, please cite:
      
! M. Simoes, C. Braithwaite, A. Makaya, E. Martinez-Paneda. Modelling  
! fatigue crack growth in Shape Memory Alloys. Fatigue & Fracture of 
! Engineering Materials & Structures 45: 1243-1257 (2022)
      
! Z. Khalil, A.Y. Elghazouli, E. Martinez-Paneda. A generalised phase 
! field model for fatigue crack growth in elastic-plastic solids with an
! efficient monolithic solver. Computer Methods in Applied Mechanics and 
! Engineering 388, 114286 (2022)
      
! Emilio Martinez-Paneda (e.martinez-paneda@imperial.ac.uk)
! Imperial College London

      module kvisual
      implicit none
      real*8 UserVar(4,11,70000) !CPE4 or CPE8R
      integer nelem
      save
      end module
      
      subroutine uelmat(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     1 props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,
     2 kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf,
     3 lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njpro,period,
     4 materiallib)

      use kvisual
      include 'aba_param.inc' !implicit real(a-h o-z)

      dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),svars(*),
     1 energy(*),coords(mcrd,nnode),u(ndofel),du(mlvarx,*),v(ndofel),
     2 a(ndofel),time(2),params(*),jdltyp(mdload,*),adlmag(mdload,*),
     3 ddlmag(mdload,*),predef(2,npredf,nnode),lflags(*),jprops(*)

      parameter(ndim=2,ntens=4,ninpt=4,nsvint=11,ndof=3) !CPE4 or CPE8R
      
      dimension wght(ninpt),dN(nnode,1),dNdz(ndim,nnode),dstran(ntens),
     1 dNdx(ndim,nnode),b(ntens,nnode*ndim),ddsdde(ntens,ntens),
     2 stress(ntens),stran(ntens),statevLocal(nsvint),coords_ip(3),
     3 predef_loc(npredf),dpredef_loc(npredf),xx1(3,3),xx1Old(3,3),
     4 olds(ntens)
      
!     preliminaries
      if(nsvars .lt. ninpt*nsvint) then
        write(7,*)'Increase the number of SDVs to', ninpt*nsvint
        call xit
      endif
      
!     find number of elements          
      if (kinc.eq.1) then
       if (jelem.eq.1) then
        nelem=jelem
       else
        if (jelem.gt.nelem) nelem=jelem
       endif 
      endif       
        
!     initialising
      if (lflags(3).eq.4) then
       do i=1,ntens
        do j=1,ntens
         ddsdde(i,j)=0.d0
        end do
        ddsdde(i,j) = 1.d0
       enddo
       do i=1,ndofel
        do j=1,ndofel
         amatrx(i,j)=0.d0
        end do
        amatrx(i,i)=1.d0
       end do
        
      else  
          
      do k1=1,ndofel
       rhs(k1,1)=0.d0
      end do
      amatrx=0.d0
      wght=1.d0
      
!     reading parameters
      xl=props(1)
      Gc=props(2)
      kflagS=props(3) ! Solution flag, 0 - monolithic, 1 - staggered
      xk=1.d-07

      do kintk=1,ninpt
!     evaluate shape functions and derivatives
       call kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)
       call kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
       call kbmatrix(dNdx,ntens,nnode,ndim,b)
       dvol=wght(kintk)*djac
       
       coords_ip=0.d0
       do k1=1,nnode
        do k2=1,mcrd
         coords_ip(k2)=coords_ip(k2)+dN(k1,1)*coords(k2,k1)
        end do
       end do
       
!     predefined fields
       if(npredf.gt.0) then
        call tempfv(kintk,ninpt,nnode,ndim,dN,predef,npredf,predef_loc,
     * dpredef_loc)
       end if
       
!     calculate phase field and incremental strains from nodal values
       phi=0.d0
       do inod=1,nnode
        phi=phi+dN(inod,1)*u(ndim*nnode+inod)
       end do
       call straininc(ntens,ndim,nnode,mlvarx,dNdx,du,dstran,u,xx1,
     * xx1Old)
       
!     recover and assign state variables       
       call kstatevar(kintk,nsvint,svars,statevLocal,1)
       stress=statevLocal(1:ntens)
       olds=stress
       stran=statevLocal((ntens+1):(2*ntens))
       phin=statevLocal(2*ntens+1)
       Hn=statevLocal(2*ntens+2)
       sed=statevLocal(2*ntens+3)
       sedn=sed

!     call umat to obtain stresses and constitutive matrix 
       dvmat=0.d0
       celent=1.d0
       call material_lib_mech(materiallib,stress,ddsdde,stran,dstran,
     *       kintk,djac,dvmat,xx1,predef_loc,dpredef_loc,npredf,celent,
     *       coords_ip)
       stran=stran+dstran
       
       if(kflagS.eq.0.or.dtime.eq.0.d0) then
        phin=phi
       end if     
       
!     compute the total strain energy density
      do i=1,ntens 
       sed=sed+(stress(i)+olds(i))*(dstran(i))/2.d0     
      end do
       
!     enforcing Karush-Kuhn-Tucker conditions
       if(kflagS.eq.0) then
        H=sed
       else
        H=sedn   
       end if    
       if (H.lt.Hn) H=Hn
       
!     store state variables        
       statevLocal(1:ntens)=stress
       statevLocal((ntens+1):(2*ntens))=stran
       statevLocal(2*ntens+1)=phi
       statevLocal(2*ntens+2)=H
       statevLocal(2*ntens+3)=sed
       call kstatevar(kintk,nsvint,svars,statevLocal,0)
       
!     form and assemble stiffness matrix and internal force vector
       amatrx(1:ndim*nnode,1:ndim*nnode)=
     1 amatrx(1:ndim*nnode,1:ndim*nnode)+
     2 dvol*(((1.d0-phin)**2+xk)*matmul(matmul(transpose(b),ddsdde),b))
        
       rhs(1:ndim*nnode,1)=rhs(1:ndim*nnode,1)-
     1 dvol*(matmul(transpose(b),stress)*((1.d0-phin)**2+xk))
     
       amatrx((ndim*nnode+1):ndofel,(ndim*nnode+1):ndofel)=
     1 amatrx((ndim*nnode+1):ndofel,(ndim*nnode+1):ndofel)+
     1 dvol*(matmul(transpose(dNdx),dNdx)*Gc*xl+
     2 matmul(dN,transpose(dN))*(Gc/xl+2.d0*H))       
  
       rhs((ndim*nnode+1):ndofel,1)=rhs((ndim*nnode+1):ndofel,1)-dvol*
     1 (matmul(transpose(dNdx),matmul(dNdx,u((ndim*nnode+1):ndofel)))*
     2 Gc*xl+dN(:,1)*((Gc/xl+2.d0*H)*phi-2.d0*H))
                   
! output
       UserVar(kintk,1:ntens,jelem)=stress(1:ntens)*((1.d0-phin)**2+xk)
       UserVar(kintk,(ntens+1):(2*ntens+3),jelem)=
     1 statevLocal((ntens+1):(2*ntens+3))
      
      end do       ! end loop on material integration points
      end if
      
      RETURN
      END

!***********************************************************************      
      subroutine kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)
c
      include 'aba_param.inc'
c
      parameter (gaussCoord=0.577350269d0)  
      dimension dN(nnode,1),dNdz(ndim,*),coord24(2,4)
      
      data  coord24 /-1.d0, -1.d0,
     2                1.d0, -1.d0,
     3               -1.d0,  1.d0,
     4                1.d0,  1.d0/  
      
      if (ninpt.eq.4.and.nnode.eq.4.and.ndim.eq.2) then ! CPE4
          
!     determine (g,h)
       g=coord24(1,kintk)*gaussCoord
       h=coord24(2,kintk)*gaussCoord

!     shape functions 
       dN(1,1)=(1.d0-g)*(1.d0-h)/4.d0
       dN(2,1)=(1.d0+g)*(1.d0-h)/4.d0
       dN(3,1)=(1.d0+g)*(1.d0+h)/4.d0
       dN(4,1)=(1.d0-g)*(1.d0+h)/4.d0

!     derivative d(Ni)/d(g)
       dNdz(1,1)=-(1.d0-h)/4.d0
       dNdz(1,2)=(1.d0-h)/4.d0
       dNdz(1,3)=(1.d0+h)/4.d0
       dNdz(1,4)=-(1.d0+h)/4.d0

!     derivative d(Ni)/d(h)
       dNdz(2,1)=-(1.d0-g)/4.d0
       dNdz(2,2)=-(1.d0+g)/4.d0
       dNdz(2,3)=(1.d0+g)/4.d0
       dNdz(2,4)=(1.d0-g)/4.d0
       
      elseif (ninpt.eq.4.and.nnode.eq.8.and.ndim.eq.2) then ! CPE8R 
          
!     determine (g,h,r)
       g=coord24(1,kintk)*gaussCoord
       h=coord24(2,kintk)*gaussCoord

!     shape functions 
       dN(1,1)=-0.25d0*(1.d0-g)*(1.d0-h)*(1.d0+g+h)
       dN(2,1)=0.25d0*(1.d0+g)*(1.d0-h)*(g-h-1.d0)
       dN(3,1)=0.25d0*(1.d0+g)*(1.d0+h)*(g+h-1.d0)
       dN(4,1)=0.25d0*(1.d0-g)*(1.d0+h)*(h-g-1.d0)
       dN(5,1)=0.5d0*(1.d0-g*g)*(1.d0-h)
       dN(6,1)=0.5d0*(1.d0+g)*(1.d0-h*h)
       dN(7,1)=0.5d0*(1.d0-g*g)*(1.d0+h)
       dN(8,1)=0.5d0*(1.d0-g)*(1.d0-h*h)        

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

      else
       write (6,*) '***ERROR: The shape fuctions cannot be found'   
      endif    
 
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
      end do

      return
      end

!***********************************************************************
      subroutine tempfv(kintk,ninpt,nnode,ndim,dN,predef,npredf,
     1 predef_loc,dpredef_loc)
c
      include 'aba_param.inc'
c
      dimension dN(nnode,1),predef(2,npredf,nnode)
      dimension predef_loc(npredf),dpredef_loc(npredf)
c
      do k1=1,npredf
       predef_loc(k1)=0.d0
       dpredef_loc(k1)=0.d0
       do k2=1,nnode
        predef_loc(k1)=predef_loc(k1)+
     &          (predef(1,k1,k2)-predef(2,k1,k2))*dN(k2,1)
        dpredef_loc(k1)=dpredef_loc(k1)+predef(2,k1,k2)*dN(k2,1)
       end do
      end do
c
      return
      end
      
!***********************************************************************
      subroutine straininc(ntens,ndim,nnode,mlvarx,bmat,du,dstran,u,xx1,
     1 xx1Old)
c
c     Notation:  
c       dstran(i)  incremental strain component 
c       note:      i = 1   xx direction 
c                    = 2   yy direction 
c                    = 3   zz direction
c                    = 4   xy direction
c    u() - displacement
c   du() - increment of displacement in the last inc.
c   
      include 'aba_param.inc'

      dimension dstran(ntens),bmat(ndim,nnode),du(mlvarx,*),xdu(3),
     1 xx1(3,3),u(mlvarx,*),utmp(3),utmpOld(3),xx1Old(3,3)

      dstran=0.d0
      ! set xx1 to Identity matrix
      xx1=0.d0
      xx1Old=0.d0
      do k1=1,3
       xx1(k1,k1)=1.d0
       xx1Old(k1,k1)=1.d0         
      end do

c************************************
c    Compute incremental strains
c************************************
c
      do nodi=1,nnode
           
       incr_row=(nodi-1)*ndim
       do i=1,ndim
        xdu(i)=du(i+incr_row,1)
        utmp(i)=u(i+incr_row,1)
        utmpOld(i)=utmp(i)-xdu(i)
       end do

       dNidx=bmat(1,nodi)
       dNidy=bmat(2,nodi)

       dstran(1)=dstran(1)+dNidx*xdu(1)
       dstran(2)=dstran(2)+dNidy*xdu(2)
       dstran(4)=dstran(4)+dNidy*xdu(1)+dNidx*xdu(2)  

c     deformation gradient
       xx1(1,1)=xx1(1,1)+dNidx*utmp(1)
       xx1(1,2)=xx1(1,2)+dNidy*utmp(1)
       xx1(2,1)=xx1(2,1)+dNidx*utmp(2)
       xx1(2,2)=xx1(2,2)+dNidy*utmp(2)
c
       xx1Old(1,1)=xx1Old(1,1)+dNidx*utmpOld(1)
       xx1Old(1,2)=xx1Old(1,2)+dNidy*utmpOld(1)
       xx1Old(2,1)=xx1Old(2,1)+dNidx*utmpOld(2)
       xx1Old(2,2)=xx1Old(2,2)+dNidy*utmpOld(2)

      end do
c
      return
      end

!***********************************************************************
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

      ddsdde=0.d0
      noffset=noel-nelem
      statev(1:nstatv)=UserVar(npt,1:nstatv,noffset)
     
      return
      end
