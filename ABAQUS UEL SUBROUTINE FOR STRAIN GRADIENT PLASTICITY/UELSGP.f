! User element subroutine for Gudmundson (2004) strain gradient 
! plasticity model. The code is distributed under a BSD license     
      
! If using this code for research or industrial purposes, please cite:
! E. Martínez-Pañeda, V.S. Deshpande, C.F. Niordson and N.A. Fleck.
! The role of plastic strain gradients in the crack growth resistance of
! metals. Journal of the Mechanics and Physics of Solids, 126: 136-150
! (2019). doi:10.1016/j.jmps.2019.02.011       
      
! Emilio Martínez-Pañeda (mail@empaneda.com)
! University of Cambridge
      
      module kvisual
      implicit none
      real*8 sigout(4,9,10000),straout(4,9,10000),epout(4,9,10000),
     * tauEout(8,9,10000),tauDout(8,9,10000),Eout(9,10000)
      save
      end module
c*****************************************************************
      subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,
     1 nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,
     2 jelem,params,ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,  
     4 ddlmag,mdload,pnewdt,jprops,njpro,period)
c
      use kvisual
      include 'aba_param.inc' !implicit real(a-h o-z)
C
      dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),svars(*),
     1 energy(*),coords(mcrd,nnode),u(ndofel),du(mlvarx,*),v(ndofel),
     2 a(ndofel),time(2),params(*),jdltyp(mdload,*),adlmag(mdload,*),
     3 ddlmag(mdload,*),predef(2,npredf,nnode),lflags(*),jprops(*)

      parameter (ndim=2, ndof=5, ndi=3, nshr=1, nnodemax=8,
     1     ntens=4, ninpt=9, nsvint=13)
c
c      ndim  ... number of spatial dimensions
c      ndof  ... number of degrees of freedom per node
c      ndi   ... number of direct stress components
c      nshr  ... number of shear stress component
c      ntens ... total number of stress tensor components (=ndi+nshr)
c      ninpt ... number of integration points
c      nsvint... number of state variables per integration pt     
c
      dimension stiff(ndof*nnodemax,ndof*nnodemax),force(ndof*nnodemax),
     1 dN(nnodemax),dshape(ndim,nnodemax),xjaci(ndim,ndim),
     2 bmat(nnodemax*ndim),statevLocal(nsvint),stress(ntens), 
     3 ddsdde(ntens,ntens),stran(ntens),dstran(ntens),wght(ninpt)
c
      dimension xjaci_new(ndim,ndim),bmat_new(nnodemax*ndim)

! SGP part
      dimension dep(ntens),gdep(ndim*ntens),qD(ntens),tauD(ndim*ntens),
     * ep(ntens),tauE(ndim*ntens),gep(ndim*ntens),tau(ndim*ntens)
c
      data wght /0.308641975309d0, 0.493827160494d0, 0.308641975309d0,
     * 0.493827160494d0, 0.79012345679d0, 0.493827160494d0, 
     * 0.308641975309d0, 0.493827160494d0, 0.308641975309d0/

!     Preliminaries
      pnewdtLocal = pnewdt
      if(jtype .ne. 1) then
        write(7,*)'Incorrect element type'
        call xit
      endif 
      if(nsvars .lt. ninpt*nsvint) then
        write(7,*)'Increase the number of SDVs to', ninpt*nsvint
        call xit
      endif 

!     Initialize rhs and lhs
      rhs(:, 1)=0.d0
      amatrx=0.d0
      ddsdde=0.d0

!     Set up elasticity matrix
      EBULK3 = props(1)/(1.d0-2.d0*props(2))
      XK = EBULK3/3.d0
      EG2 = props(1)/(1.d0+props(2))
      EG=EG2/2.d0
      ELAM = (EBULK3-EG2)/3.d0
      
      do k1 = 1, 3
       do k2 = 1, 3
        ddsdde(K2,K1) = ELAM
       end do
       ddsdde(K1,K1) = EG2 + ELAM
      end do
      ddsdde(4,4) = EG      

!     Loop over integration points
      do kintk = 1, ninpt
          
!     Initialize constitutive variables
       dep=0.d0
       gdep=0.d0
       gep=0.d0

!     Evaluate shape functions and derivatives
       call shapefcn(kintk,ninpt,nnode,ndim,dN,dshape)

!     Form B-matrix
       djac = 1.d0
       call jacobian(jelem,ndim,nnode,coords,dshape,djac,xjaci,pnewdt)
       if (pnewdt .lt. pnewdtLocal) pnewdtLocal = pnewdt

       call bmatrix(xjaci,dshape,nnode,ndim,bmat)

!     Calculate incremental strains
       call straininc(ntens,ndof,ndim,nnode,mlvarx,bmat,du,dstran,
     * u,dep,dN,gdep,gep)
       
!     Retrieve history variables and define the constitutive response       
       call statevar(kintk,nsvint,svars,statevLocal,1)
       do k1=1,ntens
        stress(k1) = statevLocal(k1)
        stran (k1) = statevLocal(k1+ntens)
        ep(k1) = statevLocal(k1+2*ntens)
       end do
       E=statevLocal(1+3*ntens)

       stress=stress+matmul(ddsdde,(dstran-dep))

!     Dissipative contributions
       dE=sqrt(2.d0/3.d0*(dep(1)**2+dep(2)**2+dep(3)**2+dep(4)**2/2.d0)
     + +props(5)**2*(gdep(1)**2+gdep(2)**2+gdep(3)**2+gdep(4)**2
     + +gdep(5)**2+gdep(6)**2+gdep(7)**2/2.d0+gdep(8)**2/2.d0))
      
       Sy=props(3)*(1.d0+(E+dE)/(props(3)/props(1)))**props(7)
       dSfdE=props(7)*props(1)*(1.d0
     + +props(1)*(E+dE)/props(3))**(props(7)-1)      
       if((E+dE).eq.0.d0) dSfdE=0.d0
      
       call kvisco(dE,dVdE,props,dtime,E,Sigma,SE,Sy)
       dSigmadE=Sy*dVdE+dSfdE*Sigma/Sy
          
       qD=2.d0/3.d0*SE*dep
       qD(4)=qD(4)/2.d0
       tauD=props(5)**2*SE*gdep
       tauD(7:8)=tauD(7:8)/2.d0
      
!     Update
       ep=ep+dep
       stran=stran+dstran
       E=E+dE
      
!     Energetic contributions  
       tauE=EG*props(4)**2*gep
       tauE(7:8)=tauE(7:8)/2.d0
       tau=tauD+tauE

       do k1=1,ntens
        statevLocal(k1) = stress(k1)
        statevLocal(ntens+k1) = stran(k1)
        statevLocal(2*ntens+k1) = ep(k1)
       end do
       statevLocal(1+3*ntens)=E

       call statevar(kintk,nsvint,svars,statevLocal,0)        

!     Visualization      
       do k1=1,ntens
        sigout(k1,kintk,jelem)=stress(k1)
        straout(k1,kintk,jelem)=stran(k1)
        epout(k1,kintk,jelem)=ep(k1)
        tauEout(k1,kintk,jelem)=tauE(k1)
        tauEout(k1+ntens,kintk,jelem)=tauE(k1+ntens)
        tauDout(k1,kintk,jelem)=tauD(k1)
        tauDout(k1+ntens,kintk,jelem)=tauD(k1+ntens)
       enddo
       Eout(kintk,jelem)=E    
        
!     Form stiffness matrix and internal force vector

       call stiffmatrix(ntens,nnode,ndim,ndof,wght(kintk),djac,ddsdde,
     1 stress,bmat,stiff,force,tau,qD,dN,dE,dep,props,gdep,Sigma,SE,
     2 dSigmadE)    

       do k1=1, ndof*nnode
        rhs(k1, 1) = rhs(k1, 1) - force(k1) 
        do k2=1, ndof*nnode
         amatrx(k1, k2) = amatrx(k1, k2) + stiff(k1,k2)
        end do
       end do
      end do       ! end loop on material integration points
      pnewdt = pnewdtLocal

      return
      end
c*****************************************************************
      subroutine shapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)
c
      include 'aba_param.inc'
c
      parameter (gaussCoord=0.774596669241483d0)
      dimension dN(*),dNdz(ndim,*),coord28(2,9)
c     
      data  coord28 /-1.d0, -1.d0,
     2                0.d0, -1.d0,
     3                1.d0, -1.d0,
     4               -1.d0,  0.d0,
     5                0.d0,  0.d0,
     6                1.d0,  0.d0,
     7               -1.d0,  1.d0,
     8                0.d0,  1.d0,      
     9                1.d0,  1.d0/
C     
C  2D 8-nodes
c
c     determine (g,h,r)
      g = coord28(1,kintk)*gaussCoord
      h = coord28(2,kintk)*gaussCoord
c
c     shape functions
      dN(1) = -0.25d0*(1.d0-g)*(1.d0-h)*(1.d0+g+h)
      dN(2) = 0.25d0*(1.d0+g)*(1.d0-h)*(g-h-1.d0)
      dN(3) = 0.25d0*(1.d0+g)*(1.d0+h)*(g+h-1.d0)
      dN(4) = 0.25d0*(1.d0-g)*(1.d0+h)*(h-g-1.d0)
      dN(5) = 0.5d0*(1.d0-g*g)*(1.d0-h)
      dN(6) = 0.5d0*(1.d0+g)*(1.d0-h*h)
      dN(7) = 0.5d0*(1.d0-g*g)*(1.d0+h)
      dN(8) = 0.5d0*(1.d0-g)*(1.d0-h*h)      
c
c     derivative d(Ni)/d(g)
      dNdz(1,1) = 0.25d0*(1.d0-h)*(2.d0*g+h)
      dNdz(1,2) = 0.25d0*(1.d0-h)*(2.d0*g-h)
      dNdz(1,3) = 0.25d0*(1.d0+h)*(2.d0*g+h)
      dNdz(1,4) = 0.25d0*(1.d0+h)*(2.d0*g-h)
      dNdz(1,5) = -g*(1.d0-h)
      dNdz(1,6) = 0.5d0*(1.d0-h*h)
      dNdz(1,7) = -g*(1.d0+h)
      dNdz(1,8) = -0.5d0*(1.d0-h*h)
c
c     derivative d(Ni)/d(h)
      dNdz(2,1) = 0.25d0*(1.d0-g)*(g+2.d0*h)
      dNdz(2,2) = 0.25d0*(1.d0+g)*(2.d0*h-g)
      dNdz(2,3) = 0.25d0*(1.d0+g)*(2.d0*h+g)
      dNdz(2,4) = 0.25d0*(1.d0-g)*(2.d0*h-g)
      dNdz(2,5) = -0.5d0*(1.d0-g*g) 
      dNdz(2,6) = -(1.d0+g)*h 
      dNdz(2,7) = 0.5d0*(1.d0-g*g)
      dNdz(2,8) = -(1.d0-g)*h     
      return
      end
c*****************************************************************
      subroutine jacobian(jelem,ndim,nnode,coords,dshape,djac,
     1 xjaci,pnewdt)
c
c     Notation: djac - Jac determinant; xjaci - inverse of Jac matrix 
c                          
      include 'aba_param.inc'
      parameter(maxDof=2)

      dimension xjac(maxDof,maxDof),xjaci(ndim,*),coords(3,*),
     * dshape(ndim,*)

      do i = 1, ndim
        do j = 1, ndim
           xjac(i,j)  = 0.d0
           xjaci(i,j) = 0.d0
        end do
      end do
c
      do inod= 1, nnode
          do idim = 1, ndim
            do jdim = 1, ndim
              xjac(jdim,idim) = xjac(jdim,idim) + 
     1        dshape(jdim,inod)*coords(idim,inod)      
            end do
          end do 
      end do

      djac = xjac(1,1)*xjac(2,2) - xjac(1,2)*xjac(2,1)
        if (djac .gt. 0.d0) then
          ! jacobian is positive - o.k.
          xjaci(1,1) =  xjac(2,2)/djac
          xjaci(2,2) =  xjac(1,1)/djac
          xjaci(1,2) = -xjac(1,2)/djac
          xjaci(2,1) = -xjac(2,1)/djac
        else
          ! negative or zero jacobian
          write(7,*)'WARNING: element',jelem,'has neg. Jacobian'
          pnewdt = 0.25d0
        endif
      return
      end
c*****************************************************************
      subroutine bmatrix(xjaci,dshape,nnode,ndim,bmat)
c
c     Notation: bmat(i) ....dN1/dx, dN1/dy, dN2/dx, dN2/dy...      
      include 'aba_param.inc'
      dimension bmat(*),dshape(ndim,*),xjaci(ndim,*)

      do i = 1, nnode*ndim
          bmat(i) = 0.d0
      end do

      do inod = 1, nnode
       do ider = 1, ndim
        do idim = 1, ndim
         irow = idim + (inod - 1)*ndim
         bmat(irow)=bmat(irow)+xjaci(idim,ider)*dshape(ider,inod) 
        end do
       end do
      end do 

      return
      end

c*****************************************************************
      subroutine straininc(ntens,ndof,ndim,nnode,
     1 mlvarx,bmat,du,dstran,u,dep,dN,gdep,gep)
c
c     Notation:  
c       dstran(i)  incremental strain component 
c       note:      i = 1   xx direction 
c                    = 2   yy direction 
c                    = 3   zz direction
c                    = 4   xy direction
c    u() - displacement / plastic strain components
c   du() - increment of displacement / etc. in the last inc.
      include 'aba_param.inc'

      dimension dstran(ntens),bmat(ndim,*),du(mlvarx,*),u(mlvarx,*),
     1 dN(nnode),xdep(3),dep(ntens),gdep(8),xdu(2),gep(8),xep(3)  

      dstran = 0.d0
c
c************************************
c    Compute incremental quantities
c************************************
c
      do nodi = 1, nnode
           
       incr_row = (nodi - 1)*ndof

       do i = 1, ndof
        if (i.le.2) then  
        xdu(i)= du(i + incr_row,1)
        else    
        xdep(i-2)=du(i + incr_row,1)
        xep(i-2)=u(i + incr_row,1)
        endif
       end do

       dNidx = bmat(1,nodi)
       dNidy = bmat(2,nodi)

       ! Strains
       dstran(1) = dstran(1) + dNidx*xdu(1)
       dstran(2) = dstran(2) + dNidy*xdu(2)
       dstran(4) = dstran(4) + dNidy*xdu(1) + dNidx*xdu(2)
         
       ! Plastic strains
       dep(1)=dep(1)+dN(nodi)*xdep(1)
       dep(2)=dep(2)+dN(nodi)*xdep(2)
       dep(3)=dep(3)-dN(nodi)*xdep(1)-dN(nodi)*xdep(2)
       dep(4)=dep(4)+dN(nodi)*xdep(3)
         
       ! Plastic strain gradient
       gdep(1)=gdep(1)+dNidx*xdep(1)
       gdep(2)=gdep(2)+dNidy*xdep(1)
       gdep(3)=gdep(3)+dNidx*xdep(2)
       gdep(4)=gdep(4)+dNidy*xdep(2)
       gdep(5)=gdep(5)-dNidx*xdep(1)-dNidx*xdep(2)
       gdep(6)=gdep(6)-dNidy*xdep(1)-dNidy*xdep(2)
       gdep(7)=gdep(7)+dNidx*xdep(3)
       gdep(8)=gdep(8)+dNidy*xdep(3)
         
       ! Total plastic strain gradient
       gep(1)=gep(1)+dNidx*xep(1)
       gep(2)=gep(2)+dNidy*xep(1)
       gep(3)=gep(3)+dNidx*xep(2)
       gep(4)=gep(4)+dNidy*xep(2)
       gep(5)=gep(5)-dNidx*xep(1)-dNidx*xep(2)
       gep(6)=gep(6)-dNidy*xep(1)-dNidy*xep(2)
       gep(7)=gep(7)+dNidx*xep(3)
       gep(8)=gep(8)+dNidy*xep(3)

      end do
c
      return
      end

c*****************************************************************
      subroutine statevar(npt,nsvint,statev,statev_ip,icopy)
c
c     Transfer data to/from element-level state variable array from/to
c     material-point level state variable array.
c
      include 'aba_param.inc'

      dimension statev(*),statev_ip(*)

      isvinc= (npt-1)*nsvint     ! integration point increment

      if (icopy .eq. 1) then
c
c     Prepare arrays for entry into umat
c
        do i = 1, nsvint
          statev_ip(i)=statev(i+isvinc)
        end do
c
      else
c
c     Update element state variables upon return from umat
c
        do i = 1, nsvint
          statev(i+isvinc)=statev_ip(i)
        end do
      end if

      return
      end

c*****************************************************************
      subroutine kvisco(dE,dVdE,props,dtime,E,Sigma,SE,Sy)    
c
c     Viscoplastic function
c
      include 'aba_param.inc'

      dimension props(*)
      
      parameter (varpi=0.01d0)
      
      xm=props(8)
      e0=props(6)

      if (nint(props(9)).eq.1) then 
      !Fuentes-Alonso and Martinez-Paneda (2019)
       e0c=e0*(1.d0/(varpi*xm))**(1/(xm-1))
       if ((dE*xm/(e0c*dtime)).le.1.d0) then
        Sigma=Sy*dE/(varpi*e0*dtime)
        dVdE=1.d0/(varpi*e0*dtime)
        SE=Sy/(varpi*e0*dtime)
       else
        Sigma=Sy*((dE/dtime-(1.d0/xm-1.d0)*e0c)/e0)**xm 
        dVdE=xm/(e0*dtime)*((dE/dtime-(1.d0/xm-1.d0)*e0c)/e0)**(xm-1.d0)
        SE=Sigma/dE
       endif      
          
      elseif (nint(props(9)).eq.2) then !Martinez-Paneda et al.(IJSS,2016)
       if (dE/(e0*dtime).le.varpi**(1.d0/xm)) then
        Sigma=Sy*dE/(e0*dtime)*varpi**(1.d0-1.d0/xm)
        dVdE=1.d0/(e0*dtime)*varpi**(1.d0-1.d0/xm)
        SE=Sy/(e0*dtime)*varpi**(1.d0-1.d0/xm)
       else
        Sigma=Sy*(dE/(e0*dtime))**xm
        dVdE=xm/(e0*dtime)*(dE/(e0*dtime))**(xm-1.d0)
        SE=Sigma/dE
       endif 
        
      else !Panteghini and Bardella(CMAME,2016)
       if (dE/(e0*dtime).le.1.d0) then
        Sigma=Sy*dE/(2.d0*e0*dtime)
        dVdE=1.d0/(2.d0*e0*dtime)
        SE=Sy/(2.d0*e0*dtime)
       else
        Sigma=Sy*(1.d0-e0*dtime/(2.d0*dE))
        dVdE=e0*dtime/(2.d0*dE**2)
        SE=Sigma/dE
       endif       
      endif
      
      return
      end

c*****************************************************************
      subroutine stiffmatrix(ntens,nnode,ndim,ndof,weight,djac,ddsdde,
     1 stress,bmat,stiff,force,tau,qD,dN,dE,dep,props,gdep,Sigma,SE,
     2 dSigmadE)    
c
c     Stiffness matrix and internal force contributions at 
c     material integration point
c
      include 'aba_param.inc'

      parameter(maxDof=5)

      dimension stiff(ndof*nnode,*),stiff_p(maxDof,maxDof),force(*),
     * force_p(maxDof),stress(ntens),bmat(ndim,*),ddsdde(ntens,ntens),
     * tau(ntens*ndim),dN(nnode),qD(ntens),dep(ntens),props(*),
     * gdep(ntens*ndim),Me(ntens,ntens),Mg(ntens*ndim,ntens*ndim) 
      dimension dEdep(1,ntens),dEgdep(1,ntens*ndim),qDdep(ntens,ntens),
     * qDgdep(ntens,ntens*ndim),tauDdep(ntens*ndim,ntens),
     * temp0(ntens,1),tauDgdep(ntens*ndim,ntens*ndim),
     * tauEgep(ntens*ndim,ntens*ndim)
      dimension Neit(3,4),Nej(4,3),Hej(8,3),temp1(4,3),temp2(8,3),
     * Heit(3,8),temp3(3,8),temp5(8,1)
	 
      double precision Me, Mg, Neit, Nej	 

      Me = 0.d0
      Mg = 0.d0
      Neit = 0.d0
      Nej = 0.d0
      Hej = 0.d0
      Heit = 0.d0

      do i=1,ntens 
       Me(i,i)=1.d0
       Mg(i,i)=1.d0
      enddo
      Me(4,4)=0.5d0
      Mg(5,5)=1.d0
      Mg(6,6)=1.d0
      Mg(7,7)=0.5d0
      Mg(8,8)=0.5d0
      tauEgep=0.5d0*props(1)/(1.d0+props(2))*props(4)**2*Mg

      do i = 1, ndof*nnode
         force(i) = 0.d0
        do j = 1, ndof*nnode
          stiff(j,i) = 0.d0
        end do
      end do

      if (dE.gt.0.d0) then   
       dEdep(1,:)=2.d0/(3.d0*dE)*matmul(Me,dep)  ! (64)
       dEgdep(1,:)=props(5)**2/dE*matmul(Mg,gdep) ! (65)
       temp=1.d0/dE*dSigmadE-Sigma/dE**2   
      else
       dSigmadE=0.d0
       dEdep(1,:)=0.d0
       dEgdep(1,:)=0.d0
       temp=0.d0
      endif
      
      temp0(:,1)=2.d0/3.d0*matmul(Me,dep)*temp
      qDdep=matmul(temp0,dEdep)+2.d0/3.d0*Me*SE ! (67)
      qDgdep=matmul(temp0,dEgdep)  ! (68)
      temp5(:,1)=props(5)**2*matmul(Mg,gdep)*temp
      tauDdep=matmul(temp5,dEdep) ! (73) 
      tauDgdep=matmul(temp5,dEgdep)+props(5)**2*Mg*SE  ! (74)
      
      dvol= weight*djac
      
      do nodj = 1, nnode
       incr_col = (nodj - 1)*ndof
       dNjdx = bmat(1,nodj)
       dNjdy = bmat(2,nodj)
       Nej(1,1)=dN(nodj)
       Nej(2,2)=dN(nodj)
       Nej(3,1)=-dN(nodj)
       Nej(3,2)=-dN(nodj)
       Nej(4,3)=dN(nodj)
       Hej(1,1)=dNjdx
       Hej(2,1)=dNjdy
       Hej(3,2)=dNjdx
       Hej(4,2)=dNjdy
       Hej(5,1)=-dNjdx
       Hej(5,2)=-dNjdx
       Hej(6,1)=-dNjdy
       Hej(6,2)=-dNjdy
       Hej(7,3)=dNjdx
       Hej(8,3)=dNjdy       
       
       force_p(1)=dNjdx*stress(1)+dNjdy*stress(4)  
       force_p(2)=dNjdy*stress(2)+dNjdx*stress(4)  
       force_p(3:5)=matmul(transpose(Nej),(qD-stress))
     * +matmul(transpose(Hej),tau)   
           
          do jdof = 1, ndof
            jcol = jdof + incr_col
            force(jcol) = force(jcol) + force_p(jdof)*dvol
          end do

          do nodi = 1, nnode
           incr_row = (nodi -1)*ndof
           dNidx = bmat(1,nodi)
           dNidy = bmat(2,nodi)
           stiff_p(1,1)=dNidx*ddsdde(1,1)*dNjdx+dNidy*ddsdde(4,4)*dNjdy
     & + dNidx*ddsdde(1,4)*dNjdy + dNidy*ddsdde(4,1)*dNjdx           
           stiff_p(1,2)=dNidx*ddsdde(1,2)*dNjdy+dNidy*ddsdde(4,4)*dNjdx
     & + dNidx*ddsdde(1,4)*dNjdx + dNidy*ddsdde(4,2)*dNjdy  
           stiff_p(2,1)=dNidy*ddsdde(2,1)*dNjdx+dNidx*ddsdde(4,4)*dNjdy
     & + dNidy*ddsdde(2,4)*dNjdy + dNidx*ddsdde(4,1)*dNjdx 
           stiff_p(2,2)=dNidy*ddsdde(2,2)*dNjdy+dNidx*ddsdde(4,4)*dNjdx
     & + dNidy*ddsdde(2,4)*dNjdx + dNidx*ddsdde(4,2)*dNjdy 
            
           stiff_p(1,3)=-dN(nodj)*(ddsdde(1,1)*dNidx+ddsdde(4,1)*dNidy
     & -ddsdde(1,3)*dNidx-ddsdde(4,3)*dNidy)               
           stiff_p(1,4)=-dN(nodj)*(ddsdde(1,2)*dNidx+ddsdde(4,2)*dNidy
     & -ddsdde(1,3)*dNidx-ddsdde(4,3)*dNidy)               
           stiff_p(1,5)=-dN(nodj)*(ddsdde(1,4)*dNidx+ddsdde(4,4)*dNidy)
           stiff_p(2,3)=-dN(nodj)*(ddsdde(2,1)*dNidy+ddsdde(4,1)*dNidx
     & -ddsdde(2,3)*dNidy-ddsdde(4,3)*dNidx)                
           stiff_p(2,4)=-dN(nodj)*(ddsdde(2,2)*dNidy+ddsdde(4,2)*dNidx
     & -ddsdde(2,3)*dNidy+ddsdde(4,3)*dNidx)              
           stiff_p(2,5)=-dN(nodj)*(ddsdde(2,4)*dNidy+ddsdde(4,4)*dNidx)
           
           stiff_p(3,1) = -dN(nodi)*(dNjdx*(ddsdde(1,1)-ddsdde(3,1))
     & +dNjdy*(ddsdde(1,4)-ddsdde(3,4)))                 
           stiff_p(3,2) = -dN(nodi)*(dNjdy*(ddsdde(1,2)-ddsdde(3,2))
     & +dNjdx*(ddsdde(1,4)-ddsdde(3,4)))                  
           stiff_p(4,1) = -dN(nodi)*(dNjdx*(ddsdde(2,1)-ddsdde(3,1))
     & +dNjdy*(ddsdde(2,4)-ddsdde(3,4)))                  
           stiff_p(4,2) = -dN(nodi)*(dNjdy*(ddsdde(2,2)-ddsdde(3,2))
     & +dNjdx*(ddsdde(2,4)-ddsdde(3,4)))                 
           stiff_p(5,1) =-dN(nodi)*(ddsdde(4,1)*dNjdx+ddsdde(4,4)*dNjdy)
           stiff_p(5,2) =-dN(nodi)*(ddsdde(4,2)*dNjdy+ddsdde(4,4)*dNjdx)
           
           Neit(1,1)=dN(nodi)
           Neit(1,3)=-dN(nodi)
           Neit(2,2)=dN(nodi)
           Neit(2,3)=-dN(nodi)
           Neit(3,4)=dN(nodi)
           Heit(1,1)=dNidx
           Heit(1,2)=dNidy
           Heit(2,3)=dNidx
           Heit(2,4)=dNidy
           Heit(1,5)=-dNidx
           Heit(2,5)=-dNidx
           Heit(1,6)=-dNidy
           Heit(2,6)=-dNidy
           Heit(3,7)=dNidx
           Heit(3,8)=dNidy        
           
           temp1=matmul((qDdep+ddsdde),Nej)+matmul(qDgdep,Hej)
           temp2=matmul(tauDdep,Nej)+matmul(tauDgdep,Hej)
           temp3=matmul(Heit,tauEgep)
           stiff_p(3:5,3:5)=matmul(Neit,temp1)+matmul(Heit,temp2)
     *      +matmul(temp3,Hej)
           
            do jdof = 1, ndof
              icol = jdof + incr_col
              do idof = 1, ndof
                irow = idof + incr_row
                stiff(irow,icol) = stiff(irow,icol) + 
     &          stiff_p(idof,jdof)*dvol
              end do
            end do
          end do
      end do
c
      return
      end

      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
      
      use kvisual
      include 'aba_param.inc'

      character*80 cmname

      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1 ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),
     2 time(2),predef(1),dpred(1),props(nprops),coords(3),drot(3,3),
     3 dfgrd0(3,3),dfgrd1(3,3)
      
      kelemshift=nint(props(1))
      ddsdde=0.0d0
      
      do k1=1,ntens
       Statev(k1)=sigout(k1,npt,noel-kelemshift)
       Statev(k1+ntens)=straout(k1,npt,noel-kelemshift)
       Statev(k1+2*ntens)=epout(k1,npt,noel-kelemshift)
       Statev(k1+3*ntens)=tauEout(k1,npt,noel-kelemshift)
       Statev(k1+4*ntens)=tauEout(k1+ntens,npt,noel-kelemshift)
       Statev(k1+5*ntens)=tauDout(k1,npt,noel-kelemshift)
       Statev(k1+6*ntens)=tauDout(k1+ntens,npt,noel-kelemshift)       
      enddo
       Statev(1+7*ntens)=Eout(npt,noel-kelemshift)      
      return
      end 