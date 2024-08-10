! HETVAL subroutines for defining diffuse interfaces
! The code is distributed under a BSD license     
      
! Adria Quintanas-Corominas (a.quintanas-corominas@imperial.ac.uk)
! Emilio Martinez-Paneda (e.martinez-paneda@imperial.ac.uk)
! Imperial College London
c***********************************************************************
      subroutine hetval(cmname,temp,time,dtime,statev,flux,predef,dpred)

      include 'aba_param.inc'
      
      character*80 cmname
      
      dimension temp(2),statev(*),predef(*),time(2),flux(2),dpred(*)

      phi=temp(1)
      xl=statev(1)
      dw=2.d0*phi
      ddw=2.d0
      flux(1)=-dw/(2.0d0*xl**2)
      flux(2)=-ddw/(2.0d0*xl**2)
      statev(2)=phi
      
      end subroutine hetval
