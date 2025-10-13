!DIR$ FREEFORM
!
! ** FIELDS
!  1: xf
!  2: xp
!  3: xb
!  4: xm
!  5: xa
!  6: Gsize [mum] (mean)
!  7: Hardness
!
! ** STATE VARIABLES
!  1: xf
!  2: xp
!  3: xb
!  4: xm
!  5: xa
!  6: Gsize [mum] (mean)
!  7: Gsize [mum]
!  8: nucf
!  9: nucp
! 10: nucb
! 11: dT_dt at 700 degrees
! 12: temp_max

! ==============================================================================
! MODULE PARAMETERS
! Parameters need to be specified in set_parameters.f90 with subroutine set_parameters

module parameters

  implicit none

  type parameterset

    ! ** Parameters to set

    ! General
    real*8 :: tolerance = 1e-4

    ! Material names
    character*80, dimension(5) :: matnames = ''

    ! Composition data
    real*8, dimension(5) :: C = 0.0, &
                            Mn = 0.0, &
                            Si = 0.0, &
                            Ni = 0.0, &
                            Cr = 0.0, &
                            Mo = 0.0, &
                            Co = 0.0, &
                            V = 0.0, &
                            W = 0.0, &
                            As = 0.0, &
                            Al = 0.0, &
                            P = 0.0, &
                            Ti = 0.0, &
                            Cu = 0.0

    ! Grain growth parameters
    real*8 :: Gsize_min = 1.0, &
              Gsize_max = 1000.0

    ! Phase-dependent thermal properties
    character*256 :: propfile = ''
    real*8, allocatable ::  data_sheat(:,:), &
                            data_cond(:,:), &
                            data_rho(:,:)

    ! Heat source parameters
    real*8 :: heat_a = 0.0, &
              heat_b = 0.0, &
              heat_cf = 0.0, &
              heat_cr = 0.0, &
              heat_Q = 0.0, &
              heat_t = 0.0     ! Set equal to step time
    logical, dimension(10) :: heat_step = .false.
    real*8, dimension(2,10) :: heat_center = 0.0


    ! ** Derived parameters composition (can also be specified)
    real*8, dimension(5) :: Ae1 = 0.d0, &
                            Ae3 = 0.d0, &
                            Bs = 0.d0,  &
                            Ms = 0.d0,  &
                            Me = 0.d0,  &
                            Ac3 = 0.d0, &
                            fcomp_f, &
                            fcomp_p, &
                            fcomp_b, &
                            xf_eq, &
                            Hvfp, Hvfplog, &
                            Hvb, Hvblog, &
                            Hvm, Hvmlog


    !** Derived parameters heat source 
    real*8 :: heat_v, &
              heat_tau, &
              heat_ff, &
              heat_fr

  end type

  type(parameterset) :: param

end module

! ==============================================================================
! SUBROUTINE UEXTERNALDB
! Used to set and calculate parameters at the beginning of the analysis

subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)

    use parameters

    implicit none 

    ! Input
    integer, intent(in) :: lop, lrestart, kstep, kinc
    real*8, intent(in) :: time(2), dtime

    ! Local
    real*8 :: Ae1_grange, Ae1_andrews, Ae1_eldis, &
      Ae3_grange, Ae3_andrews, Ae3_eldis, &
      c_gamma, c_gamma_grange, c_gamma_andrews, c_gamma_eldis
    integer :: imat

    ! Pre-calculate at start of analysis
    if( lop == 0 ) then

      include 'set_parameters_thermal.f90'

      do imat = 1, 5 

        if( param%Ae3(imat) == 0.0 ) then
          Ae3_grange = (1570. - 323.*(param%C(imat)) - 25.*param%Mn(imat) + 80.*param%Si(imat) - 32.*param%Ni(imat) - 3.*param%Cr(imat) - 30.)*5./9.
          Ae3_andrews = 910. - 203.*sqrt((param%C(imat))) + 44.7*param%Si(imat) - 15.2*param%Ni(imat) + 31.5*param%Mo(imat) + 104.*param%V(imat) &
              + 13.1*param%W(imat) - 30.*param%Mn(imat) + 11.*param%Cr(imat) + 20.*param%Cu(imat) - 700.*param%P(imat) - 400.*param%Al(imat) - 120.*param%As(imat)- 400.*param%Ti(imat)
          Ae3_eldis = 871.0 - 254.4*sqrt((param%C(imat))) + 51.7*param%Si(imat) - 14.2*param%Ni(imat)
          param%Ae3(imat) = (Ae3_grange + Ae3_andrews + Ae3_eldis) / 3.0
        end if
        if( param%Ae1(imat) == 0.0 ) then
          Ae1_grange = (1333. - 25.*param%Mn(imat) + 40.*param%Si(imat) - 26.*param%Ni(imat) + 42.*param%Cr(imat) - 32.)*5./9.
          Ae1_andrews = 723. - 16.9*param%Ni(imat) + 29.1*param%Si(imat) + 6.38*param%W(imat) - 10.7*param%Mn(imat) + 16.9*param%Cr(imat) + 290*param%As(imat)
          Ae1_eldis = 712. - 17.8*param%Mn(imat) + 20.1*param%Si(imat) - 19.1*param%Ni(imat) + 11.9*param%Cr(imat) + 9.8*param%Mo(imat)
          param%Ae1(imat) = (Ae1_grange + Ae1_andrews + Ae1_eldis) / 3.0
        end if
        if( param%Bs(imat) == 0.0 ) then
          param%Bs(imat) = 637. - 58.*(param%C(imat)) - 35.*param%Mn(imat) - 15.*param%Ni(imat) - 34.*param%Cr(imat) - 41.*param%Mo(imat)
        end if
        if( param%Ms(imat) == 0.0 ) then
          param%Ms(imat) = 539. - 423.*(param%C(imat)) - 30.4*param%Mn(imat) - 17.7*param%Ni(imat) - 12.1*param%Cr(imat) - 7.5*param%Mo(imat) + 10.*param%Co(imat) - 7.5*param%Si(imat)
        end if
        if( param%Me(imat) == 0.0 ) then
          param%Me(imat) = 80.0
        end if
        if( param%Ac3(imat) == 0.0 ) then
          param%Ac3(imat) = param%Ae3(imat)
        end if

        if( param%fcomp_f(imat) == 0.0 ) then
          param%fcomp_f(imat) = exp(1.0 + 6.31*(param%C(imat)) + 1.78*param%Mn(imat) + 0.31*param%Si(imat) + 1.12*param%Ni(imat) + 2.7*param%Cr(imat) + 4.06*param%Mo(imat))
        end if
        if( param%fcomp_p(imat) == 0.0 ) then
          param%fcomp_p(imat) = exp(-4.25 + 4.12*(param%C(imat)) + 4.36*param%Mn(imat) + 0.44*param%Si(imat) + 1.71*param%Ni(imat) + 3.33*param%Cr(imat) + 5.19*param%Mo(imat)**0.5)
        end if
        if( param%fcomp_b(imat) == 0.0 ) then
          param%fcomp_b(imat) = exp(-10.23 + 10.18*(param%C(imat)) + 0.85*param%Mn(imat) + 0.55*param%Ni(imat) + 0.90*param%Cr(imat) + 0.36*param%Mo(imat))
        end if

        if( param%xf_eq(imat) == 0.0 ) then
          c_gamma_grange = (-param%Ae1(imat)*9./5. + 1570. - 25.*param%Mn(imat) + 80.*param%Si(imat) - 32.*param%Ni(imat) - 3.*param%Cr(imat) - 30.)/323.
          c_gamma_andrews = ((-param%Ae1(imat) + 910. - 203.*sqrt((param%C(imat))) + 44.7*param%Si(imat) - 15.2*param%Ni(imat) + 31.5*param%Mo(imat) + 104.*param%V(imat) &
              + 13.1*param%W(imat) - 30.*param%Mn(imat) + 11.*param%Cr(imat) + 20.*param%Cu(imat) - 700.*param%P(imat) - 400.*param%Al(imat) - 120.*param%As(imat)- 400.*param%Ti(imat)) &
              / 203.)**2
          c_gamma_eldis = ((-param%Ae1(imat) + 871.0 + 51.7*param%Si(imat) - 14.2*param%Ni(imat))/254.4)**2
          c_gamma = (c_gamma_grange + c_gamma_andrews + c_gamma_eldis) / 3.0
          param%xf_eq(imat) = (c_gamma - (param%C(imat))) / c_gamma
        end if

        param%Hvfp(imat) = 42. + 223.*(param%C(imat)) + 53.*param%Si(imat) + 30.*param%Mn(imat) + 12.6*param%Ni(imat) + 7.*param%Cr(imat) + 19.*param%Mo(imat)
        param%Hvfplog(imat) = 10. - 19.*param%Si(imat) + 4.*param%Ni(imat) + 8.*param%Cr(imat) + 130.*param%V(imat)
        param%Hvb(imat) = -323. + 185.*(param%C(imat)) + 330.*param%Si(imat) + 153.*param%Mn(imat) + 65.*param%Ni(imat) + 144.*param%Cr(imat) + 191.*param%Mo(imat)
        param%Hvblog(imat) = 89. + 53.*(param%C(imat)) - 55.*param%Si(imat) - 22.*param%Mn(imat) - 10.*param%Ni(imat) - 20.*param%Cr(imat) - 33.*param%Mo(imat)
        param%Hvm(imat) = 127. + 949.*(param%C(imat)) + 27.*param%Si(imat) + 11.*param%Mn(imat) + 8.*param%Ni(imat) + 16.*param%Cr(imat)
        param%Hvmlog(imat) = 21.

        write(6,*) ''
        write(6,*) '===== Thermodynamic values ====='
        write(6,*) 'Material', imat
        write(6,*) 'Ae3', param%Ae3(imat)
        write(6,*) 'Ae1', param%Ae1(imat)
        write(6,*) 'Bs ', param%Bs(imat)
        write(6,*) 'Ms ', param%Ms(imat)
        write(6,*) 'xf_eq', param%xf_eq(imat)
        write(6,*) '================================'
        write(6,*) ''

      end do

      ! Speed and delay heat source
      ! midpoint starts and ends at 5% of maximum power
      param%heat_v = (param%heat_cf + param%heat_cr) / param%heat_t
      param%heat_tau = param%heat_cf*param%heat_t / (param%heat_cf + param%heat_cr)
      param%heat_ff = 2.0*param%heat_cf / (param%heat_cf + param%heat_cr)
      param%heat_fr = 2.0 - param%heat_fr

      ! Read properties
      call k_read_props()

    end if

end subroutine

! ==============================================================================
! SUBROUTINE USDFLD
! Main function of the phase transformation calculation

subroutine usdfld(field,statev,pnewdt,direct,t,celent, &
  time,dtime,cmname,orname,nfield,nstatv,noel,npt,layer, &
  kspt,kstep,kinc,ndi,nshr,coord,jmac,jmatyp,matlayo,laccfla)

  use parameters
  implicit none

  ! Input
  character*80, intent(in) :: cmname, orname
  real*8, intent(in) :: direct(3,3), t(3,3), celent, time(2), dtime, coord(*)
  integer, intent(in) :: nfield, nstatv, noel, npt, layer, kspt, kstep, kinc, ndi, nshr, &
     jmac(*), jmatyp(8), matlayo, laccfla

  ! In/output 
  real*8, intent(inout) :: statev(nstatv), pnewdt

  ! Output
  real*8, intent(out) :: field(nfield)

  ! Local
  character*3  :: flgray(15)
  integer :: imat, jarray(15), jrcd
  real*8 :: array(15)

  ! Get material index
  do imat = 1, 5
    if( param%matnames(imat) == cmname ) then 
      exit
    end if 
  end do
  
  if( imat == 6 ) then
    write(6,*) ''
    write(6,*) 'MATERIAL NAME'
    write(6,*) cmname
    write(6,*) 'NOT FOUND'
    write(6,*) ''
    call xit
  end if

  ! Initial step, initialize state variables
  if( sum(statev(1:5)) < 0.9 ) then

    ! Get temperature
    call getvrm('TEMP',array,jarray,flgray,jrcd,jmac,jmatyp,matlayo,laccfla)
  
    if( array(1) > param%Ac3(imat) ) then
      statev(1:4) = 0.d0
      statev(5) = 1.d0
      statev(8:10) = 0.d0
    else
      statev(1) = min(max(min(field(1), param%xf_eq(imat)), 0.0), 1.0)
      if( statev(1) > param%tolerance * param%xf_eq(imat) ) then
        statev(8) = 1.1
      else
        statev(8) = 0.0
      end if

      statev(2) = min(max(min(field(2), 1.0-param%xf_eq(imat)), 0.0), 1.0)
      if( statev(2) > param%tolerance * (1.0-param%xf_eq(imat)) ) then
        statev(9) = 1.1
      else
        statev(9) = 0.0
      end if

      statev(3) = min(max(field(3), 0.0), 1.0)
      if( statev(3) > param%tolerance ) then
        statev(10) = 1.1
      else
        statev(10) = 0.0
      end if

      statev(4) = min(max(field(4), 0.0), 1.0)
      statev(5) = max(1.0 - sum(statev(1:4)), 0.0)
    end if

    ! Initial grain size
    statev(6:7) = max(field(6), param%Gsize_min)

  end if

  field(1:6) = statev(1:6)
  
  ! Calculate hardness
  if( statev(11) > 1e-12 ) then
    field(7) = (field(1)+field(2))*(param%Hvfp(imat)+param%Hvfplog(imat)*log10(statev(11)*3600)) &
      + field(3)*(param%Hvb(imat)+param%Hvblog(imat)*log10(statev(11)*3600)) &
      + field(4)*(param%Hvm(imat)+param%Hvmlog(imat)*log10(statev(11)*3600))
  end if

end subroutine

! ==============================================================================
! SUBROUTINE UMATHT

subroutine umatht(u,dudt,dudg,flux,dfdt,dfdg,&
        statev,temp,dtemp,dtemdx,time,dtime,predef,dpred,&
        cmname,ntgrd,nstatv,props,nprops,coords,pnewdt,&
        noel,npt,layer,kspt,kstep,kinc,vold,co,lakonl,konl,&
        ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,mi)
  
  use parameters
  implicit none
  
  character*8 lakonl
  character*80 cmname
  
  integer :: ntgrd,nstatv,nprops,noel,npt,layer,kspt,kstep,kinc,&
    konl(20),ipompc(*),nodempc(3,*),nmpc,ikmpc(*),ilmpc(*),mi(*)
  
  real*8 u,dudt,dudg(ntgrd),flux(ntgrd),dfdt(ntgrd),&
    statev(nstatv),pnewdt,temp,dtemp,dtemdx(ntgrd),time(2),dtime,&
    predef(*),dpred(*),props(nprops),coords(3),dfdg(ntgrd,ntgrd),&
    vold(0:mi(2),*),co(3,*),coefmpc(*)
  
  integer :: i, imat
  real*8 :: condp, sheatp, rhop
  real*8 :: cond, sheat, rho
  real*8 :: temp_prev, dtemp_prev, dt_tr, temp_tr, fun_tc, &
            Ac1, Ac3, Ae1, Ae3, Bs, Ms, Me, &
            corr, xf, xp, xb, xm, xa, dxf, dxp, dxb, dxm, dxa, nucf, nucp, nucb, &
            Gsize, dGsize, Gastm, Gsize_mean

  ! Get material index
  do imat = 1, 5
    if( param%matnames(imat) == cmname ) then 
      exit
    end if 
  end do
  
  if( imat == 6 ) then
    write(6,*) ''
    write(6,*) 'MATERIAL NAME'
    write(6,*) cmname
    write(6,*) 'NOT FOUND'
    write(6,*) ''
    call xit
  end if

  ! Get statevars
  xf = statev(1)
  xp = statev(2)
  xb = statev(3)
  xm = statev(4)
  xa = statev(5)
  Gsize_mean = statev(6)
  Gsize = statev(7)
  nucf = statev(8)
  nucp = statev(9)
  nucb = statev(10)

  ! Temperature increase
  temp_prev = temp
  temp = temp + dtemp

  ! Set inital values
  dxf = 0.0
  dxp = 0.0
  dxb = 0.0
  dxm = 0.0
  dxa = 0.0
  corr = 1.0  

  ! === Grain growth ===
  if( dtemp >= 0.0 .and. temp > param%Ae1(imat) ) then

    ! First austenization step > reset initial grain size
    if( temp_prev <= param%Ae1(imat) ) then
      Gsize = xa*Gsize + (1.0 - xa)*param%Gsize_min
    end if

    call k_split_inc_heat(dt_tr, temp_tr, dtime, temp, temp_prev, param%Ae1(imat))
    call k_increment_grain_size(dGsize, Gsize, temp_tr, dt_tr)
    Gsize = Gsize + dGsize

  else if( temp <= param%Ae3(imat) .and. temp_prev > param%Ae3(imat) ) then

    ! Average grain size for partial austenization
    Gsize_mean = xa*Gsize + (1.0 - xa)*Gsize_mean

  end if

  ! === Heating stage ===
  if( dtemp > 0.0 .and. xa < 1.0-param%tolerance ) then
    ! Austenite
    if( temp > param%Ae3(imat) ) then
      call k_split_inc_heat(dt_tr, temp_tr, dtime, temp, temp_prev, param%Ae3(imat))
      call k_inc_aus(dxa, xa, temp_tr, param%Ae1(imat), param%Ae3(imat), dt_tr)
      dxf = - dxa * xf / (1.0 - xa)
      dxp = - dxa * xp / (1.0 - xa)
      dxb = - dxa * xb / (1.0 - xa)
      dxm = - dxa * xm / (1.0 - xa)
    end if

  ! === Cooling stage ===
  else if( dtemp < 0.0 .and. xa > 0.0+param%tolerance ) then
    ! ASTM Grain size number
    Gastm = 2.0 * log(254.0/Gsize) / log(2.0) + 1.0

    ! Ferrite
    if( (temp > param%Bs(imat) .and. temp <= param%Ae3(imat)) .or. (temp_prev > param%Bs(imat) .and. temp_prev <= param%Ae3(imat)) ) then
      call k_split_inc_cool(dt_tr, temp_tr, dtime, temp, temp_prev, param%Ae3(imat), param%Bs(imat))
      fun_tc = (param%Ae3(imat) - temp_tr)**3 * exp(-1.384e4/(temp_tr+273.15)) * 2.0**(0.41*Gastm) / param%fcomp_f(imat)
      call k_inc_diff(dxf, nucf, xf/param%xf_eq(imat), fun_tc, temp_tr, dt_tr)
      dxf = dxf * param%xf_eq(imat)
    end if

    ! Pearlite
    if( (temp > param%Bs(imat) .and. temp <= param%Ae1(imat)) .or. (temp_prev > param%Bs(imat) .and. temp_prev <= param%Ae1(imat)) ) then
      call k_split_inc_cool(dt_tr, temp_tr, dtime, temp, temp_prev, param%Ae1(imat), param%Bs(imat))
      fun_tc = (param%Ae1(imat) - temp_tr)**3 * exp(-1.384e4/(temp_tr+273.15)) * 2.0**(0.32*Gastm) / param%fcomp_p(imat)
      call k_inc_diff(dxp, nucp, xp/(1.0-param%xf_eq(imat)), fun_tc, temp_tr, dt_tr)
      dxp = dxp * (1.0-param%xf_eq(imat))
    end if

    ! Bainite
    if( (temp > param%Ms(imat) .and. temp <= param%Bs(imat)) .or. (temp_prev > param%Ms(imat) .and. temp_prev <= param%Bs(imat)) ) then
      call k_split_inc_cool(dt_tr, temp_tr, dtime, temp, temp_prev, param%Bs(imat), param%Ms(imat))
      fun_tc = (param%Bs(imat) - temp_tr)**2 * exp(-1.384e4/(temp_tr+273.15)) * 2.0**(0.29*Gastm) / param%fcomp_b(imat)
      call k_inc_diff(dxb, nucb, xb, fun_tc, temp_tr, dt_tr)
    end if

    ! Martensite
    if( (temp > param%Me(imat) .and. temp <= param%Ms(imat)) .or. (temp_prev > param%Me(imat) .and. temp_prev <= param%Ms(imat)) ) then
      call k_split_inc_cool(dt_tr, temp_tr, dtime, temp, temp_prev, param%Ms(imat), param%Me(imat))
      call inc_mar(dxm, xm, xa, temp_tr, param%Ms(imat))
    end if

    dxa = -(dxf + dxp + dxb + dxm)
  end if

  ! Correct for over/under shooting
  if( xa + dxa > 1.0 ) then
    corr = (1.0 - xa) / dxa
  elseif( xa + dxa < 0.0 ) then
    corr = - xa / dxa
  end if

  ! Set new phase fractions
  xf = xf + corr * dxf
  xp = xp + corr * dxp
  xb = xb + corr * dxb
  xm = xm + corr * dxm
  xa = xa + corr * dxa

  ! Phases have to re-nucleate if almost 0
  if( xf < param%tolerance * param%xf_eq(imat) .and. nucf >= 1.0 ) then
    nucf = 0.0
  end if
  if( xp < param%tolerance * (1.0-param%xf_eq(imat)) .and. nucp >= 1.0 ) then
    nucp = 0.0
  end if
  if( xb < param%tolerance .and. nucb >= 1.0 ) then
    nucb = 0.0
  end if

  ! update field and statevars
  statev(1) = min(max(xf, 0.0), 1.0)
  statev(2) = min(max(xp, 0.0), 1.0)
  statev(3) = min(max(xb, 0.0), 1.0)
  statev(4) = min(max(xm, 0.0), 1.0)
  statev(5) = min(max(xa, 0.0), 1.0)
  statev(6) = Gsize_mean
  statev(7) = Gsize
  statev(8) = nucf
  statev(9) = nucp
  statev(10) = nucb

  ! Set CR at 700 degrees
  if( temp <= 700 .and. temp_prev > 700 ) then
    statev(11) = abs(dtemp) / dtime
  else if( temp > 700 ) then
    statev(11) = 0.0
  end if

  ! Maximum temperature
  statev(12) = max(statev(12), temp)

  ! Material behavior
  call k_interp_prop(cond, temp+dtemp, statev(5), param%data_cond, size(param%data_cond,1), .false.)
  call k_interp_prop(sheat, temp+dtemp, statev(5), param%data_sheat, size(param%data_sheat,1), .false.)
  call k_interp_prop(rho, temp+dtemp, statev(5), param%data_rho, size(param%data_rho,1), .true.)
  
  dudt = sheat * rho
  u = u + dudt * dtemp

  do i = 1, ntgrd
    flux(i) = -cond * dtemdx(i)
    dfdg(i,i) = -cond
  end do 

end subroutine

! ==============================================================================
! SUBROUTINE k_inc_aus
! Calculates the incremenet in austenite phase fraction upon heating with
! Leblond-Devaux model

subroutine k_inc_aus(dx, x, temp, Ae1, Ae3, dtime)

  implicit none

  ! Input/output
  real*8, intent(in) :: x, temp, Ae1, Ae3, dtime
  real*8, intent(out) :: dx

  ! Local
  real*8 :: x_eq, tau, r, dr_dx

  if( temp <= Ae3 ) then
    x_eq = (temp-Ae1)/(Ae3-Ae1)
    tau = 1.0 - x_eq*0.8
  else
    x_eq = 1.0
    tau = 0.05
  end if

  dx = dtime * (x_eq - x) / (tau + dtime)

end subroutine

! ==============================================================================
! SUBROUTINE k_inc_diff
! Calculates the incremenet interf phase fraction for diffuse transofmrations
! upon heating with the KV-model of Li et al. 1998

subroutine k_inc_diff(dx, nuc, x_in, fun_tc, temp, dtime)

  use parameters

  implicit none

  ! Input/output
  real*8, intent(in) :: x_in, fun_tc, temp, dtime
  real*8, intent(inout) :: nuc
  real*8, intent(out) :: dx

  ! Local
  real*8 :: x, dt, dx_nuc, dx_dt, r, x1, x2, x3, &
    r1, r2, r3, cr, cs, ct, cp, cq
  integer :: iter

  dt = dtime
  dx_nuc = 0.0
  dx = 0.0

  ! Not yet nucleated
  if( nuc < 1.0 ) then

    dx_dt = fun_tc / 0.10434035495809084
    nuc = nuc + dx_dt * dt

    ! Remaining time after nucleation
    if( nuc >= 1.0 ) then
      dx_nuc = 0.01
      dt = (nuc - 1.0) / dx_dt
    end if

  end if

  ! Nucleated
  if( nuc >= 1.0 ) then

    x1 = x_in + dx_nuc
    x3 = 1.0
    x2 = 0.5 * (x1 + x3)
 
    r1 = x1 - x_in - dx_nuc - dt * fun_tc * x1**(0.4*(1.0-x1)) * (1.0-x1)**(0.4*x1)
    r2 = x2 - x_in - dx_nuc - dt * fun_tc * x2**(0.4*(1.0-x2)) * (1.0-x2)**(0.4*x2)
    r3 = x3 - x_in - dx_nuc - dt * fun_tc * x3**(0.4*(1.0-x3)) * (1.0-x3)**(0.4*x3)

    iter = 0

    ! Brent root finding algorithm
    do while(.true.)
      iter = iter + 1

      cr = r2 / r3
      cs = r2 / r1
      ct = r1 / r3  
      cp = cs*(ct*(cr-ct)*(x3-x2) - (1.0-cr)*(x2-x1))
      cq = (ct-1.0)*(cr-1.0)*(cs-1.0)

      x = x2 + cp / cq

      ! Bisection method if outside bounds
      if( .not.((x > x1 .and. x < x2 .and. r1*r2 < 0.0) &
          .or. (x > x2 .and. x < x3 .and. r2*r3 < 0.0)) ) then 
        if( (r1 < 0.0 .and. r2 > 0.0) .or. (r1 > 0.0 .and. r2 < 0.0) ) then 
          x = 0.5 * (x1 + x2)
        else 
          x = 0.5 * (x2 + x3)
        end if
      end if

      ! Residual
      r = x - x_in - dx_nuc - dt * fun_tc * x**(0.4*(1.0-x)) * (1.0-x)**(0.4*x)

      if( abs(r) < param%tolerance ) then
        exit
      elseif( iter > 100 ) then
        write(6,*) ''
      	write(6,*) 'DIFFUSE INCREMENT DID NOT CONVERGE'
      	write(6,*) ''
      	write(6,*) r, x_in, x
      	call xit
      end if

      ! Update points 
      if( x > x1 .and. x < x2 ) then 
        x3 = x2 
        r3 = r2
      else 
        x1 = x2 
        r1 = r2
      end if 

      x2 = x 
      r2 = r
      
    end do

    dx = x - x_in

  end if

end subroutine

! ==============================================================================
! SUBROUTINE INC_MAR
! Calculates the incremenet in martensite upon heating using the
! Koistenen-Marburger equation

subroutine inc_mar(dx, x, xa, temp, Ms)

  implicit none

  real*8, intent(out) :: dx
  real*8, intent(in) :: x, xa, temp, Ms

  dx = (xa + x) * (1.0 - exp(-1.1e-2*(Ms-temp))) - x

end subroutine

! ==============================================================================
! SUBROUTINE k_split_inc_cool
! Splits the part of the time of the increment in the transformation temperature
! upon cooling

subroutine k_split_inc_cool(dt_tr, T_tr, dtime, T, T_prev, Ts, Te)

  implicit none

  real*8, intent(out) :: dt_tr, T_tr
  real*8, intent(in) :: dtime, T, T_prev, Ts, Te

  if( T_prev > Ts ) then
    dt_tr = (Ts-T)/(T_prev-T)*dtime
    T_tr = T
  elseif( T < Te ) then
    dt_tr = (T_prev-Te)/(T_prev-T)*dtime
    T_tr = Te
  else
    dt_tr = dtime
    T_tr = T
  end if

end subroutine

! ==============================================================================
! SUBROUTINE k_split_inc_heat
! Splits the part of the time of the increment in the transformation temperature
! upon heating

subroutine k_split_inc_heat(dt_tr, T_tr, dtime, T, T_prev, Ts)

  implicit none

  real*8, intent(out) :: dt_tr, T_tr
  real*8, intent(in) :: dtime, T, T_prev, Ts

  if( T_prev < Ts ) then
    dt_tr = (T-Ts)/(T-T_prev)*dtime
  else
    dt_tr = dtime
  end if
  T_tr = T

end subroutine

! ==============================================================================
! SUBROUTINE k_increment_grain_size
! Calculates  increment the austenite grain size based on
! Pous-Romero et al. 2013

subroutine k_increment_grain_size(dGsize, Gsize, temp, dt)

  use parameters
  implicit none

  real*8, intent(out) :: dGsize
  real*8, intent(in) :: Gsize, temp, dt

  dGsize = dt * 2.4e8 * exp(-22853.0/(temp+273.15)) * (1.0/Gsize - 1.0/param%Gsize_max)
  dGsize = dt * 10.0e8 * exp(-22853.0/(temp+273.15)) * (1.0/Gsize - 1.0/param%Gsize_max)

end subroutine

! ==============================================================================
! SUBROUTINE INTERP_PROPS
! Calculates phase and temperature dependent properties
! NOTE: Somehow asume-size arrays do not work with Abaqus?????

subroutine k_interp_prop(prop, temp, xa, table, ndata, inverse)

  implicit none

  ! Input
  real*8, intent(out) :: prop
  real*8, intent(in) :: temp, xa
  real*8, dimension(ndata,3), intent(in) :: table
  integer, intent(in) :: ndata
  logical, intent(in) :: inverse

  ! Local
  integer :: i
  real*8 :: phase_prop(2), phase_dpropdt(2)

  if( temp <= table(1,3) ) then
    phase_prop = table(1,1:2)

  elseif( temp >= table(ndata,3) ) then
    phase_prop = table(ndata,1:2)
    
  else
    do i = 1, ndata - 1
      if ( temp >= table(i,3) .and. temp < table(i+1,3)) then
        phase_prop = table(i,1:2) &
          + (temp - table(i,3))/(table(i+1,3) - table(i,3)) &
          * (table(i+1,1:2) - table(i,1:2))
        exit
      end if
    end do
  end if

  if( inverse ) then
    prop = 1.d0/((1.d0 - xa)/phase_prop(1) + xa/phase_prop(2))
  else
    prop = (1.d0 - xa)*phase_prop(1) + xa*phase_prop(2)
  end if

end subroutine

! ==============================================================================
! SUBROUTINE READ_PROPS
! Read phase dependent properties from file

subroutine k_read_props

  use parameters 
  implicit none

  character*256 :: header
  integer :: npoints, ipoint, iprop
  logical :: exist

  inquire(file=trim(param%propfile), exist=exist)
  if( .not. exist ) then 
    write(6,*) 'INPUT FILE NOT FOUND'
    call xit 
  end if

  iprop = 0
  open(unit=101, file=trim(param%propfile))
  do while( iprop /= 3 )
    read(101,*) header
    read(101,*) npoints

    if( header == "*conductivity" ) then
      allocate(param%data_cond(npoints,3))
      do ipoint = 1, npoints
         read(101,*) param%data_cond(ipoint,:)
      end do
    elseif( header == "*specificheat" ) then
      allocate(param%data_sheat(npoints,3))
      do ipoint = 1, npoints
         read(101,*) param%data_sheat(ipoint,:)
      end do
    elseif( header == "*density" ) then
      allocate(param%data_rho(npoints,3))
      do ipoint = 1, npoints
         read(101,*) param%data_rho(ipoint,:)
      end do
    else 
      write(6,*) 'INPUTFILE HEADER IS WRONG'
      write(6,*) header
      call xit 
    end if

    iprop = iprop + 1
  end do
  close(unit=101)
end subroutine k_read_props
