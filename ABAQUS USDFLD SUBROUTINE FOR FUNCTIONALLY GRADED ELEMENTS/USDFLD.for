************************************************************************
**  USDFLD FOR FGMs        *********************************************
************************************************************************
** Emilio Martinez-Paneda
** University of Cambridge    
** mail@empaneda.com
** If using this code for research or industrial purposes, please cite:
** Emilio Martinez-Paneda, On the finite element implementation of 
** functionally graded materials, Materials 12(2): 287 (2019)

      subroutine usdfld(field,statev,pnewdt,direct,t,celent,
     1 time,dtime,cmname,orname,nfield,nstatv,noel,npt,layer,
     2 kspt,kstep,kinc,ndi,nshr,coord,jmac,jmatyp,matlayo,laccfla)
c
      include 'aba_param.inc'
c
      character*80 cmname,orname
      character*3  flgray(15)
      dimension field(nfield),statev(nstatv),direct(3,3),
     1 t(3,3),time(2)
      dimension array(15),jarray(15),jmac(*),jmatyp(*),coord(*)
c    dimension intv(1),realv(1)

      x=coord(1)
      y=coord(2)
      
      field(1)=dexp(2.07944154168*x)

      return
      end
      