! User element subroutine for quadratic plane strain cohesive elements
! Includes viscous regularization (Gao and Bower, 2004) and both
! cyclic and monotonic damage - see details in the attached document.
! The code is distributed under a BSD license     
      
! If using this code for research or industrial purposes, please cite:
! S. del Busto, C. Betegón, E. Martínez-Pañeda. A cohesive zone 
! framework for environmentally assisted fatigue. Engineering Fracture
! Mechanics 185: 210-226 (2017). doi:10.1016/j.engfracmech.2017.05.021     
      
! Emilio Martínez-Pañeda (mail@empaneda.com)
! Technical University of Denmark

      subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     1 props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,
     2 kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf,
     3 lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njpro,period)
     
      include 'aba_param.inc' !implicit real(a-h o-z)
      
      dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),svars(*),
     1 energy(*),coords(mcrd,nnode),u(ndofel),du(mlvarx,*),v(ndofel),
     2 a(ndofel),time(2),params(*),jdltyp(mdload,*),adlmag(mdload,*),
     3 ddlmag(mdload,*),predef(2,npredf,nnode),lflags(*),jprops(*)
	 
      dimension Sc(ndofel,ndofel),Fc(ndofel,nrhs),V_l(ndofel),
     1 T(mcrd,nrhs),T_d(mcrd,mcrd),U_l(ndofel),R(mcrd,mcrd),
     2 Bc(mcrd,ndofel),Bct(ndofel,mcrd),ShapeN(nnode),del(mcrd),
     3 tmp(ndofel,mcrd),delD(2),DelT(3*mcrd),DispSepM(ndofel,ndofel),
     4 Rb(ndofel,ndofel),ShapeM(mcrd,ndofel),DelTi(3*mcrd),SN(3) 
	 
!     Step-1: Read input data & Initialize
      kflag=props(6)
      th=1.d0
      n_GP=props(8)
      rhs(:,1)=0.d0
      amatrx=0.d0
      T=0.d0
      Rb=0.d0
      ShapeM=0.d0
      DispSepM=0.d0
      do i=1,6 
       DispSepM(i,i)=1.d0
      enddo

!     Step-2: Obtain the local separation from the global nodal disp.
      call kCoordTrans(R,el_length,coords,u,ndofel,nnode,mcrd,kflag)
      
      do i=1,mcrd
       do j=1,mcrd
        do k=0,nnode-1
         Rb(i+2*k,j+2*k)=R(i,j)
        enddo
       enddo
      enddo
      
      do i=0,nnode-1
       U_l(1+i*mcrd)=R(1,1)*U(1+i*mcrd)+R(1,2)*U(2+i*mcrd)
       U_l(2+i*mcrd)=R(2,1)*U(1+i*mcrd)+R(2,2)*U(2+i*mcrd)
       V_l(1+i*mcrd)=R(1,1)*V(1+i*mcrd)+R(1,2)*V(2+i*mcrd)
       V_l(2+i*mcrd)=R(2,1)*V(1+i*mcrd)+R(2,2)*V(2+i*mcrd)		 
      enddo

      DelT=matmul(DispSepM,U_l)
      DelTi=matmul(DispSepM,V_l)

      do i=1,n_GP
       D0=svars(i)
       call kgauss(n_GP,GP,GP_w,i)
       SN(1)=(-GP/2.d0)*(1.d0-GP) 
       SN(2)=(1.d0-GP)*(1.d0+GP)  
       SN(3)=(GP/2.d0)*(1.d0+GP) 

       if (D0.ge.1.d0) then
        T=0.d0
        T_d=0.d0
       else               
        if((kflag.eq.0) .or. (kflag.eq.1)) then         
         del(1)=SN(1)*DelT(1)+SN(2)*DelT(3)+SN(3)*DelT(5)
         delD(1)=SN(1)*DelTi(1)+SN(2)*DelTi(3)+SN(3)*DelTi(5)
        elseif(kflag.eq.2) then         
         del(1)=0.d0
         delD(1)=0.d0
        endif
        del(2)=SN(1)*DelT(2)+SN(2)*DelT(4)+SN(3)*DelT(6)
        delD(2)=SN(1)*DelTi(2)+SN(2)*DelTi(4)+SN(3)*DelTi(6)
            
        delm0 = svars(i+12)
        Tm0 = svars(i+24) 
        delt0 = svars(i+36)
        deltT = svars(i+48)
        DeltaS=4*props(2)
        Cf=props(7)

! Determine load case
        if (del(2).ge.delm0) then ! uploading
         call kSeplaw(props,del,delD,T,T_d,dtime,delm0,Tm0,D0,0)
         svars(i+12)=del(2)
         svars(i+24)=T(2,1)
        elseif ((del(2).lt.0.d0).and.(D0.ge.1.d0)) then ! contact
         call kSeplaw(props,del,delD,T,T_d,dtime,delm0,Tm0,D0,1)
        elseif ((del(2).lt.0.d0).and.(D0.lt.1.d0)) then ! compression
         call kSeplaw(props,del,delD,T,T_d,dtime,delm0,Tm0,D0,2)   
        else ! unloading/reloading
         call kSeplaw(props,del,delD,T,T_d,dtime,delm0,Tm0,D0,3)    
        endif

! Damage accumulation
        deltT=deltT+abs(del(2)-delt0)
        if (deltT.gt.props(2)) then
         Dc=abs(del(2)-delt0)/DeltaS*(T(2,1)/(props(1)*(1.d0-D0))-Cf)
         if (Dc.lt.0.d0) then
          Dc=0.d0
         endif
        else
         Dc=0.d0  
        endif
        if(delm0.gt.props(2)) then
         Dm=(del(2)-delm0)/(4*props(2))
        else
         Dm=0.d0
        endif
        D=D0+max(Dc,Dm)
        if (D.ge.1.d0) D=1.d0

        svars(i)=D
        svars(i+36)=del(2)
        svars(i+48)=deltT
       endif
!
!     Step-4: Obtain the nodal forces and the stiffess matrix
       ShapeM(1,1)=SN(1)
       ShapeM(2,2)=SN(1)
       ShapeM(1,3)=SN(2)
       ShapeM(2,4)=SN(2)
       ShapeM(1,5)=SN(3)
       ShapeM(2,6)=SN(3)         
         
       Bc=matmul(ShapeM,Rb)
       Bct=transpose(Bc)
       tmp=matmul(Bct,T_d)
       Sc=matmul(tmp,Bc)
       Fc=matmul(Bct,T)
       thick=0.5d0*el_length*GP_w*th
       call kMatrixScalar(amatrx,Sc,thick,ndofel,ndofel)
       call kMatrixScalar(rhs,-Fc,thick,ndofel,nrhs)
      enddo
 
      RETURN
      END   
! End of the main subrutine

! Complementary subroutines:
!
! Subroutine kCoordsTrans to obtrain R, the two dimensional
! transformation matrix and el_length, the element size.
      subroutine kCoordTrans(R,el_length,coords,U,ndofel,
     & nnode,mcrd,kflag)
      include 'aba_param.inc'
      dimension R(mcrd,mcrd),COORDS(mcrd,nnode),U(ndofel)
      dimension Co_de(mcrd,nnode),Co_de_m(2,2)
      
! Computation of the nodal coordinates in the deformed configuration.
      if((kflag.eq.0).or.(kflag.eq.2)) then      
       do i=1,mcrd
        do j=1,nnode
         Co_de(i,j)=COORDS(i,j)+U(2*(j-1)+i)
        end do
       end do
      elseif(kflag.eq.1) then
       Co_de=COORDS
      endif        

! Calculate of the directional cosine & the transformation matrix
      d_x=Co_de(1,3)-Co_de(1,1)
      d_y=Co_de(2,3)-Co_de(2,1)
      el_length=(d_x**2+d_y**2)**0.5d0
      cos_a=d_x/el_length
      sin_a=d_y/el_length
      R(1,1)=cos_a
      R(1,2)=sin_a
      R(2,1)=-sin_a
      R(2,2)=cos_a   
      return
      end
! End of kCoordsTransform subroutine

! Subroutine kMatrix_PlusScalar to multiply a matrix
! with a scalar number
      subroutine kMatrixScalar(A,B,c,n,m)
      include 'aba_param.inc'
      dimension A(n,m),B(n,m)
      do i=1,n
       do j=1,m
        A(i,j)=A(i,j)+c*B(i,j)
       enddo
      enddo
      return
      end
! End of subroutine kMatrix_PlusScalar

! Subroutine to assign Gauss weights and coordinates
      subroutine kgauss(n_GP,GP,GP_w,i)
      include 'aba_param.inc'
       
      if (n_GP.eq.3) then
       if(i.eq.1) then
        GP=-0.774596669241483d0
        GP_w=0.55555555556d0
       elseif(i.eq.2) then
        GP=0.d0
        GP_w=0.88888888889d0     
       elseif(i.eq.3) then
        GP=0.774596669241483d0
        GP_w=0.55555555556d0              
       endif
      elseif (n_GP.eq.6) then   
       if(i.eq.1) then
        GP=-0.932469514203152d0
        GP_w=0.1713244923791709d0
       elseif(i.eq.2) then
        GP=-0.6612093864662646d0
        GP_w=0.3607615730481379d0  
       elseif(i.eq.3) then
        GP=-0.2386191860831968d0
        GP_w=0.4679139345726913d0
       elseif(i.eq.4) then
        GP=0.2386191860831968d0
        GP_w=0.4679139345726913d0    
       elseif(i.eq.5) then
        GP=0.6612093864662646d0
        GP_w=0.3607615730481379d0
       elseif(i.eq.6) then
        GP=0.932469514203152d0
        GP_w=0.1713244923791709d0   
       endif
      elseif (n_GP.eq.12) then   
       if(i.eq.1) then
        GP=-0.981560634246732d0
        GP_w=0.04717533638647547d0
       elseif(i.eq.2) then
        GP=-0.904117256370452d0
        GP_w=0.1069393259953637d0
       elseif(i.eq.3) then
        GP=-0.7699026741943177d0
        GP_w=0.1600783285433586d0
       elseif(i.eq.4) then
        GP=-0.5873179542866143d0
        GP_w=0.2031674267230672d0
       elseif(i.eq.5) then
        GP=-0.3678314989981804d0
        GP_w=0.2334925365383534d0
       elseif(i.eq.6) then
        GP=-0.12523340851114688d0
        GP_w=0.2491470458134027d0
       elseif(i.eq.7) then
        GP=0.12523340851114688d0
        GP_w=0.2491470458134027d0
       elseif(i.eq.8) then
        GP=0.3678314989981804d0
        GP_w=0.2334925365383534d0
       elseif(i.eq.9) then
        GP=0.5873179542866143d0
        GP_w=0.2031674267230672d0   
       elseif(i.eq.10) then
        GP=0.7699026741943177d0
        GP_w=0.1600783285433586d0
       elseif(i.eq.11) then
        GP=0.904117256370452d0
        GP_w=0.1069393259953637d0    
       elseif(i.eq.12) then
        GP=0.981560634246732d0
        GP_w=0.04717533638647547d0          
       endif
      endif
      RETURN
      END          

! Subroutine kSeplaw accounting for the traction-separation law
      subroutine kSeplaw(props,del,delD,T,T_d,dtime,delm0,Tm0,D0,kflag)
!
      include 'aba_param.inc'
!
      dimension props(9),del(2),delD(2),T(2,1),T_d(2,2)
! Currently coded for Xu-Needleman constitutive law, see:
! X.P. Xu and A. Needleman, Modelling Simul. Mater. Sci. Eng. 1, 111-132

! del(1) & T(1,1): tangential separation and traction
! del(2) & T(2,1): normal separation and traction
      
      dn=props(2)
      dt=props(3)
      q=props(4)
      r=props(5)      
      Smax=props(1)*(1.d0-D0)
      if (D0.ge.1.0) Smax=0.d0
      sepwrk=dexp(1.d0)*Smax*dn
      alph=10.d0
      
      if (kflag.eq.0) then !uploading
       T(2,1)=(sepwrk/dn)*dexp(-del(2)/dn)*(del(2)/dn)
       T_d(2,2)=(sepwrk/(dn*dn))*dexp(-del(2)/dn)*(1.d0-del(2)/dn) 
      elseif (kflag.eq.1) then ! contact
       T(2,1)=alph*props(1)*dexp(1.d0)*dexp(-del(2)/dn)*del(2)/dn
       T_d(2,2)=alph*props(1)*dexp(1.d0)/dn*dexp(-del(2)/dn)*
     * (1-del(2)/dn)
      elseif (kflag.eq.2) then ! compression
       T_d(2,2)=11.d0*Smax*exp(1.d0)/dn*dexp(-del(2)/dn)*(1-del(2)/dn)
       T(2,1)=Smax*dexp(1.d0)*del(2)/dn*dexp(-del(2)/dn)+Tm0-Smax*
     * dexp(1.d0)*delm0/dn+10.d0*props(1)*dexp(1.d0)*del(2)/dn*
     * dexp(-del(2)/dn)
      else !unloading/reloading
       T_d(2,2)=Smax*dexp(1.d0)/dn
       T(2,1)=Tm0+T_d(2,2)*(del(2)-delm0)
      endif
      
      ! Mode I
      T(1,1)=0.d0
      T_d(1,1)=0.d0
      T_d(1,2)=0.d0
      T_d(2,1)=0.d0  

      ! Viscous regularization (Gao and Bower, MSMSE 2004)
      xi=props(9)
      T(2,1)=T(2,1)+xi*Smax*delD(2)/dn
      T_d(2,2)=T_d(2,2)+xi*Smax/dn/dtime   
!
      return
      end
! End of subroutine kSeplaw