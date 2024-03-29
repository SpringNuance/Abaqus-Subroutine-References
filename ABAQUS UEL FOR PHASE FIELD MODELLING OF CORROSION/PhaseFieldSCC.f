! User element subroutine for the phase field model for corrosion    
! Plane strai version (CPE4 and CPE8R elements)
! The code is distributed under a BSD license     
      
! If using this code for research or industrial purposes, please cite:
! C. Cui, R. Ma and E. Martinez-Paneda, 
! A phase field formulation for dissolution-driven stress corrosion cracking. 
! Journal of the Mechanics and Physics of Solids 147, 104254 (2021)
! doi: https://doi.org/10.1016/j.jmps.2020.104254.

      module kvisual
      implicit none
      real*8 UserVar(4,24,70000),PitD
      integer nelem
      save
      end module
      
!***********************************************************************
      subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)
      include 'aba_param.inc' 
      !implicit real(a-h o-z)
      dimension time(2)
      if (lop.eq.0) then !start of analysis
       call mutexinit(1)
      endif
      return
      end 
      
!***********************************************************************
      subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     1 props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,
     2 kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf,
     3 lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njpro,period)

      use kvisual
      include 'aba_param.inc' 
      !implicit real(a-h o-z)
      
      dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),svars(*),
     1 energy(*),coords(mcrd,nnode),u(ndofel),du(mlvarx,*),v(ndofel),
     2 a(ndofel),time(2),params(*),jdltyp(mdload,*),adlmag(mdload,*),
     3 ddlmag(mdload,*),predef(2,npredf,nnode),lflags(*),jprops(*)

      parameter(ndim=2,ntens=4,ninpt=4,nsvint=24,ndof=4)
      
      dimension wght(ninpt),dN(nnode,1),dNdz(ndim,nnode),stress(ntens),
     1 dNdx(ndim,nnode),b(ntens,nnode*ndim),ddsdde(ntens,ntens),
     2 stran(ntens),xm(nnode,nnode),BB(nnode,nnode),dstran(ntens),
     3 statevLocal(nsvint),eelas(ntens)
        
!     initialising
      do k1=1,ndofel
       rhs(k1,1)=0.d0
      end do
      amatrx=0.d0
      wght=1.d0
      
!     find number of elements          
      if (dtime.eq.0.d0) then
       if (jelem.eq.1) then   
        nelem=jelem
       else
        CALL MutexLock(1)
        if (jelem.gt.nelem) nelem=jelem 
        CALL MutexUnlock(1)
       endif 
      endif      
      
!     reading parameters
      D=props(5)
      xL0=props(6)
      a_phi=props(7)
      wh=props(8)
      AA=props(9)
      xk=props(10)
      ef=props(11)
      t0=props(12)
      cse=1.d0        ! Normalized equilibrium concentration of solid
      cle=0.036d0     ! Normalized equilibrium concentration of liquid
      xkap=1.d-7      ! well-conditioning parameter
      T=300.d0        ! Temperature (K)
      R=8314.d0       ! Gas constant   
      
      do kintk=1,ninpt
!     evaluate shape functions and derivatives
       call kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)
       call kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
       call kbmatrix(dNdx,ntens,nnode,ndim,b)
       dvol=wght(kintk)*djac                   
       
!     compmute phase field, concentration and strains from nodal values
       phi=0.d0
       dphi=0.d0
       cL=0.d0
       do inod=1,nnode
        phi=phi+dN(inod,1)*u(ndim*nnode+inod)
        dphi=dphi+dN(inod,1)*du(ndim*nnode+inod,1)
        cL=cL+dN(inod,1)*u((ndim+1)*nnode+inod)
       end do
       
       dstran=matmul(b,du(1:ndim*nnode,1))
       
!     defining suitable functions       
       hphi=-2.d0*phi**3+3.d0*phi**2
       dhphi=-6.d0*phi**2+6.d0*phi
       ddhphi=-12.d0*phi+6.d0
       gphi=(1.d0-phi)*(1.d0-phi)*phi**2
       dgphi=2.d0*phi+4.d0*phi**3-6.d0*phi**2
       ddgphi=2.d0+12.d0*phi**2-12.d0*phi
           
!     recover and assign state variables
       call kstatevar(kintk,nsvint,svars,statevLocal,1)
       stress=statevLocal(1:ntens)
       stran(1:ntens)=statevLocal((ntens+1):(2*ntens))
       eelas(1:ntens)=statevLocal((2*ntens+1):3*ntens)
       ti=statevLocal(4*ntens+7)
       ei=statevLocal(4*ntens+8)

!     update SCC length
       coordy=0.13d0
       if (phi.le.0.5d0) then
        do i=1,nnode
         coordy=coordy-dN(i,1)*coords(2,i)
        enddo
!        coordy=0.13d0-coordy
        if (coordy.gt.PitD) then
         PitD=coordy
        endif 
       endif
       if (time(2).lt.3.01d0) PitD=0.d0       
       
!     call umat to obtain stresses and constitutive matrix
       call kumat(props,ddsdde,stress,dstran,ntens,statevLocal,eelas,
     1 spd,hphi,phi,xkap,Sh0,eqplas,deqpl)
       
!     enhance the value of L 
       Sh=(hphi+xkap)*Sh0
       pls=statevLocal(4*ntens+1)
       xkm=dexp(Sh*7.12d3/(R*T))*(1.d0+pls/(props(3)/props(1)))
       if (time(2).lt.3.01d0) then
        ti=0.d0
        ei=0.d0
       else
        ei=ei+statevLocal(4*ntens+5)
        ti=ti+dtime        
       endif

       if (ei.gt.ef) then
        ti=0.d0
        ei=0.d0
       endif

       if (ti.lt.t0) then       
        xL=xL0*xkm
       else    
        xL=xL0*xkm*dexp(-xk*(ti-t0))
       endif

!     store state variables
       statevLocal(1:ntens)=stress(1:ntens)
       statevLocal((ntens+1):(2*ntens))=statevLocal((ntens+1):(2*ntens))
     1 +dstran(1:ntens)
       statevLocal(4*ntens+2)=phi
       statevLocal(4*ntens+3)=xL
       statevLocal(4*ntens+4)=Sh
       statevLocal(4*ntens+6)=cL
       statevLocal(4*ntens+7)=ti
       statevLocal(4*ntens+8)=ei
       call kstatevar(kintk,nsvint,svars,statevLocal,0) 

!     form and assemble stiffness matrix and internal force vector
       
       amatrx(1:16,1:16)=amatrx(1:16,1:16)+dvol*
     2 ((hphi+xkap)*matmul(matmul(transpose(b),ddsdde),b))
        
       rhs(1:16,1)=rhs(1:16,1)-
     1 dvol*(matmul(transpose(b),stress)*(hphi+xkap))        
       
       fc=AA*((cL-hphi*(cse-cle)-cle)**2)+wh*gphi
       dfc=-2.d0*AA*(cL-hphi*(cse-cle)-cle)*(cse-cle)*dhphi+wh*dgphi
       ddfc=wh*ddgphi-2.d0*AA*(cse-cle)*
     1 (ddhphi*(cL-hphi*(cse-cle)-cle)-dhphi*dhphi*(cse-cle))
       BB=matmul(transpose(dNdx),dNdx)
       xm=matmul(dN,transpose(dN))

       amatrx(17:24,17:24)=amatrx(17:24,17:24)+
     1 +dvol*(xm/dtime/xL+xm*ddfc+a_phi*BB) 

       rhs(17:24,1)=rhs(17:24,1)-dvol*
     1 (matmul(xm,du(17:24,1))/dtime/xL+dN(:,1)*dfc
     2 +a_phi*matmul(transpose(dNdx),matmul(dNdx,u(17:24))))
       
       ww=(cse-cle)*dhphi

       amatrx(25:32,25:32)=amatrx(25:32,25:32)+dvol*(xm/dtime+D*BB)

       rhs(25:32,1)=rhs(25:32,1)-dvol*(
     1 matmul(xm,du(25:32,1))/dtime+D*(matmul(transpose(dNdx),
     2 matmul(dNdx,u(25:32)))
     3 -ww*matmul(transpose(dNdx),matmul(dNdx,u(17:24)))))
              
! output
       UserVar(kintk,1:ntens,jelem)=stress(1:ntens)*(hphi+xkap)
       UserVar(kintk,(ntens+1):(4*ntens+8),jelem)=
     1 statevLocal((ntens+1):(4*ntens+8))
!       UserVar(kintk,4*ntens+9,jelem)=PitD
      
      end do       ! end loop on material integration points
      
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
           
          
!     determine (g,h)
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
      subroutine kumat(props,ddsdde,stress,dstran,ntens,statev,eelas,
     1 spd,hphi,phi,xkap,Sh,eqplas,deqpl)
c
c     Subroutine with the material model
c
      include 'aba_param.inc' !implicit real(a-h o-z)
      
      dimension props(*),ddsdde(ntens,ntens),stress(ntens),statev(*),
     + dstran(ntens),eelas(ntens),eplas(ntens),flow(ntens),olds(ntens),
     + oldpl(ntens) 
      
      parameter(toler=1.d-6,newton=20)  

      eelas(1:ntens)=statev((2*ntens+1):3*ntens)
      eplas(1:ntens)=statev((3*ntens+1):4*ntens)
      eqplas=statev(1+4*ntens)
      deqpl=statev(5+4*ntens)
      olds=stress
      oldpl=eplas
      
!     Initialization
      ddsdde=0.d0
      E=props(1) ! Young's modulus
      xnu=props(2) ! Poisson's ratio
      Sy=props(3) ! Yield stress
      xn=props(4) ! Strain hardening exponent    
      
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
      
!     Get yield stress from the specified Ehardening curve
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
!        if(deqpl.lt.0.d0) deqpl=-deqpl
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
      Sh=(stress(1)+stress(2)+stress(3))/3.d0
      statev((2*ntens+1):3*ntens)=eelas
      statev((3*ntens+1):4*ntens)=eplas
      statev(4*ntens+1)=eqplas 
      statev(4*ntens+5)=deqpl

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
      statev(1:nstatv-1)=UserVar(npt,1:nstatv-1,noffset)
      statev(nstatv)=PitD
      return
      end
