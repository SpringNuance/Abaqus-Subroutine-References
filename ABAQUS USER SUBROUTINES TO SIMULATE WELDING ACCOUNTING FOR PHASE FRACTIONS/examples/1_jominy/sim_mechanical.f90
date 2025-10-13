!DIR$ FREEFORM

! ** PREDEFINED FIELDS
!  1: xf
!  2: xp
!  3: xb
!  4: xm
!  5: xa
!
! ** STATE VARIABLES
!  1:     THE   (thermal strain)
!  2:     PEEQ  (eq total plastic strain; plastic+TRIP)
!  3:     PPEEQ (eq plastic strain)
!  4:     TPEEQ (eq TRIP strain)
!  5-8:   PE    (total plastic strain)
!  9:12:  PPE   (plastic strain)
!  13:16: TPE   (TRIP strain)
!  17:20: EE    (elastic strain)
!  21:    SPD   (plastic dissipation)
!  22:    Y     (yield stress, including hardening)
!  23:    Y0    (initial yield stress)
!  24:    K     (bulk modulus)
!  25:    G     (shear modulus)
!  26:    RHO   (density)
!  27:    A     (TRIP parameter)
!  28:    RHOR  (reference density)
!  29:    SYF   (base factor for sy)
!  30:    PEEQT (plastic strain total, not annealed)

! ==============================================================================
! MODULE PARAMETERS

module parameters

  implicit none

  type parameterset

    ! ** Parameters to set

    real*8 :: tolerance = 1d-4

    real*8 :: T_anneal = 5000d0
    real*8 :: n = 0d0

    character(len=256) :: propfile
    real*8, dimension(:,:), allocatable :: data_shear, &
                                           data_bulk, &
                                           data_yield, &
                                           data_trip, &
                                           data_density

    real*8 :: sy_factor = 1.d0

  end type

  type(parameterset) :: param

end module

! ==============================================================================
! MODULE MATH

module math

  implicit none

  ! ** Unit tensors in Voigt format
  real*8, dimension(4,4), parameter :: I4v = reshape([ &
    1d0, 0d0, 0d0, 0d0, &
    0d0, 1d0, 0d0, 0d0, &
    0d0, 0d0, 1d0, 0d0, &
    0d0, 0d0, 0d0, .5d0 &
  ], [4,4])

  real*8, dimension(4,4), parameter :: I2I2v = reshape([ &
    1d0, 1d0, 1d0, 0d0, &
    1d0, 1d0, 1d0, 0d0, &
    1d0, 1d0, 1d0, 0d0, &
    0d0, 0d0, 0d0, 0d0 &
  ], [4,4])

  real*8, dimension(4,4), parameter :: I4dv = I4v - I2I2v / 3.0

  ! ** Pre-calculated square-roots
  real*8, parameter :: sqrt32 = sqrt(3.d0/2.d0), &
                       sqrt23 = sqrt(2.d0/3.d0)

  contains

    ! ** FUNCTION OUTER
    function outer(v1, v2) result(v1v2)

      ! Input
      real*8, dimension(:), intent(in) :: v1, v2

      ! Output
      real*8, dimension(size(v1),size(v2)) :: v1v2

      ! Local
      integer :: i, j

      do concurrent( i=1:size(v1), j=1:size(v2) )
        v1v2(i,j) = v1(i) * v2(j)
      end do

    end function

end module math

! ==============================================================================
! SUBROUTINE UEXTERNALDB
! Used to set and calculate parameters at the beginning of the analysis

subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)

  use parameters
  implicit none

  integer :: lop, lrestart, kstep, kinc
  real*8 :: time(2), dtime

  ! Pre-calculate at start of analysis
  if( lop == 0 ) then

    include 'set_parameters_mechanical.f90'

    call read_props

  end if

end subroutine

! ==============================================================================
! SUBROUTINE UEXPAN
! Thermal expansion and phase transformation strain

subroutine uexpan(expan,dexpandt,temp,time,dtime,predef, &
       dpred,statev,cmname,nstatv,noel)

  use parameters
  implicit none

  character*80 cmname

  real*8 :: expan(*),dexpandt(*),temp(2),time(2),dtime, predef(*), &
    dpred(*),statev(nstatv)
  integer :: nstatv, noel

  ! local
  real*8 :: values(5), predefb(4), dpredb(4), eps_th, rho_a, rho_i, rho_mean, xa, dxa
  
  ! Bound phase fractions (correct extrapolation errors)
  predefb = min(max(predef(1:4), 0.d0), 1.d0)
  dpredb = dpred(1:4) + predefb - predef(1:4)
  xa = (1.d0-sum(predefb))
  dxa = -sum(dpredb)

  ! Average initial yield stress
  call interp_prop(values, temp(1), param%data_yield, 5, size(param%data_yield,1))
  statev(23) = sum(predefb * values(1:4)) + xa * values(5)

  ! Average Bulk modulus
  call interp_prop(values, temp(1), param%data_bulk, 5, size(param%data_bulk,1))
  statev(24) = sum(predefb * values(1:4)) + xa * values(5)

  ! Average Shear modulus
  call interp_prop(values, temp(1), param%data_shear, 5, size(param%data_shear,1))
  statev(25) = sum(predefb * values(1:4)) + xa * values(5)

  ! Average density
  call interp_prop(values, temp(1), param%data_density, 5, size(param%data_density,1))
  statev(26) = 1.d0 / (sum(predefb / values(1:4)) + xa / values(5))

  ! Total TRIP parameter
  if( dxa < 0.d0 ) then
    call interp_prop(values(1:4), temp(1), param%data_trip, 4, size(param%data_trip,1))
    statev(27) = max(sum(1.5d0 * values(1:4) * (2.d0-2.d0*predefb) * dpredb), 0.d0)
  else
    statev(27) = 0.d0
  end if

  ! Reference expansion in first step
  if( statev(28) == 0.d0 ) then
    if( noel == 6778 ) then 
      write(6,*) 'Element initiated'
      write(6,*) 'temp', temp  
      write(6,*) 'xa', xa
      write(6,*) 'predefb', predefb
      write(6,*) 'time', time
      write(6,*) 'rho', statev(26)
      write(6,*) 'rho phases', values(1:5)
    end if
    statev(28) = statev(26)
  end if

  ! Yield stress base factor
  if( statev(29) == 0.d0 ) then 
    statev(29) = param%sy_factor
  else
    statev(29) = min(1.0d0 + (1.d0-xa)*(param%sy_factor-1.0d0), statev(29))
  end if

  ! Current expansion
  eps_th = (statev(28) / statev(26))**(1.d0/3.d0) - 1.d0

  ! Thermal strain increment
  expan(1) = eps_th - statev(1)

  ! Update statevars
  statev(1) = eps_th

end

! ==============================================================================
! SUBROUTINE UMAT
! Thermal expansion and phase transformation strain

subroutine umat(stress,statev,ddsdde,sse,spd,scd, &
       rpl,ddsddt,drplde,drpldt, &
       stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname, &
       ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt, &
       celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

  use parameters
  use math
  implicit none

  character*80 cmname
  integer :: ntens, nstatv, ndi, nshr, nprops, noel, npt, layer, kspt, kstep, kinc
  real*8 :: sse, spd, scd, pnewdt, dtime, rpl, drpldt, temp, dtemp, celent
  real*8 :: stress(ntens),statev(nstatv), &
   ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens), &
   stran(ntens),dstran(ntens),time(2),predef(*),dpred(*), &
   props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

  real*8, dimension(4) :: eps, eps_p_tot, eps_e, eps_ed, sig_d, n
  real*8 :: eps_ev, sig_v, sig_eq_tr, A, sig_eq, lam_tot, dlam, dsigy_ddlam, &
    sigy0, sigy, bulkmod, shearmod, eps_ptp
  logical :: anneal, failed

  ! Total strain
  eps = stran + dstran

  ! Get statevars
  eps_p_tot = statev(5:8)
  lam_tot = statev(2)
  sigy0 = statev(23)*statev(29)
  bulkmod = statev(24)
  shearmod = statev(25)
  A = statev(27)

  ! Elastic strain, volumetric-deviatoric split
  eps_e = eps - eps_p_tot
  eps_ev = sum(eps_e(1:3))
  eps_ed = eps_e - eps_ev / 3.d0 * [1.d0, 1.d0, 1.d0, 0.d0]

  ! Trial stress
  sig_v = bulkmod * eps_ev
  sig_d = 2.d0 * shearmod * eps_ed * [1.d0, 1.d0, 1.d0, 0.5d0]
  sig_eq_tr = sqrt(1.5d0 * (sig_d(1)**2 + sig_d(2)**2 + sig_d(3)**2 + 2.d0*sig_d(4)**2))

  ! Initial plastic multip increment
  dlam = 0.d0
  failed = .false.

  ! Get hardened yield stress
  if( temp > param%T_anneal ) then
    anneal = .true.
    sigy = sigy0
  else 
    anneal = .false.
    call get_hardening(sigy, dsigy_ddlam, dlam, sigy0, shearmod, bulkmod, 0.d0, lam_tot, failed)
  end if
  statev(22) = sigy

  ! Only elastic deformation
  if( sig_eq_tr < sigy .and. (A == 0.d0 .or. sig_eq_tr == 0.d0)) then
    stress = sig_v * [1.d0, 1.d0, 1.d0, 0.d0] + sig_d
    ddsdde = bulkmod * I2I2v + 2.d0*shearmod * I4dv
    statev(17:20) = eps_e
    return
  end if

  !******************
  !** Plasticity part

  ! Check for yielding
  if( sig_eq_tr > sigy ) then
    call plastic_return_map(sig_eq, dlam, dsigy_ddlam, sigy0, sig_eq_tr, shearmod, bulkmod, A, lam_tot, anneal, failed)
    sigy = sig_eq
  end if

  ! Only TRIP
  if( dlam <= 0.d0 ) then
    dlam = 0.d0
    call trip_return_map(sig_eq, sig_eq_tr, shearmod, A, failed)

    if( .not. anneal ) then
      call get_hardening(sigy, dsigy_ddlam, dlam, sigy0, shearmod, bulkmod, 0.d0, lam_tot + 2.d0/3.d0*A*sig_eq, failed)
    end if
  end if

  ! Return if not converged
  if( failed ) then
    pnewdt = 0.1d0
    return
  end if

  ! Plastic inc direction
  n = sqrt32 * sig_d / sig_eq_tr

  ! Update stress
  sig_d = sqrt23 * sig_eq * n
  stress = sig_v * [1.d0, 1.d0, 1.d0, 0.d0] + sig_d

  ! Update strains
  statev(9:12) = statev(9:12) + sqrt32 * dlam * n * [1.d0, 1.d0, 1.d0, 2.d0]
  statev(13:16) = statev(13:16) + A * sig_d * [1.d0, 1.d0, 1.d0, 2.d0]
  statev(5:8) = statev(9:12) + statev(13:16)
  statev(17:20) = eps - statev(5:8)

  ! Update equivalent strains
  if( anneal ) then 
    statev(2:4) = 0.d0
    statev(21) = 0.d0
  else
    statev(3) = statev(3) + dlam
    statev(4) = statev(4) + 2.d0/3.d0 * A * sig_eq
    statev(2) = statev(3) + statev(4)
    statev(21) = statev(21) + (dlam + 2.d0/3.d0*A*sig_eq) * sig_eq
  end if
  statev(30) = statev(30) + dlam + 2.d0/3.d0*A*sig_eq

  ! Update hardened yield stress
  statev(22) = sigy

  ! Tangent
  if( dlam > 0.d0 ) then
    ddsdde = bulkmod* I2I2v &
      + 2.d0*shearmod * (1.d0 - 3.d0*shearmod*dlam/sig_eq_tr - 2.d0*shearmod*A*sig_eq/sig_eq_tr) * I4dv &
      + (6.d0*shearmod**2 * (dlam/sig_eq_tr - 1.d0/(3.d0*shearmod + dsigy_ddlam*(1.d0 + 2.d0*shearmod*A))) &
      +  4.d0*shearmod**2 * (A*sig_eq/sig_eq_tr - A*dsigy_ddlam/(3.d0*shearmod + dsigy_ddlam*(1.d0 + 2.d0*shearmod*A))) &
        ) * outer(n, n)

  ! Only TRIP tangent
  else
    ddsdde = bulkmod* I2I2v &
      + 2.d0*shearmod * (1.d0 - 2.d0*shearmod*A*sig_eq/sig_eq_tr) * I4dv &
      + 4.d0*shearmod**2 * A*sig_eq/sig_eq_tr * outer(n, n)
  end if  

end subroutine

! ==============================================================================
! SUBROUTINE PLASTIC_RETURN_MAP
! Return map algorithm for plasticity with TRIP

subroutine plastic_return_map(sig_eq, dlam, dsigy_ddlam, sigy0, sig_eq_tr, shearmod, bulkmod, A, lam_tot, anneal, failed)

  use parameters
  implicit none

  ! Output
  real*8, intent(out) :: sig_eq, dlam, dsigy_ddlam

  ! Input
  real*8, intent(in) :: sigy0, sig_eq_tr, shearmod, bulkmod, A, lam_tot
  logical, intent(in) :: anneal
  
  ! Convergence flag
  logical, intent(inout) :: failed

  ! Local
  integer :: iter
  real*8 :: r, dr_ddlam

  dlam = 0.d0
  iter = 0

  do

    if( anneal ) then 
      sig_eq = sigy0
      dsigy_ddlam = 0.d0
    else
      call get_hardening(sig_eq, dsigy_ddlam, dlam, sigy0, shearmod, bulkmod, A, lam_tot, failed)
      if( failed ) then; return; end if
    end if

    r = sig_eq_tr - 3.d0*shearmod*dlam - (1.d0 + 2.d0*shearmod*A)*sig_eq

    if( abs(r) < param%tolerance * sig_eq_tr ) then
      exit
    end if

    if( iter == 100 ) then
      write(7,*) ''
      write(7,*) '*** PLASTIC RETURN MAP DID NOT CONVERGE'
      failed = .true.
      return
    end if

    dr_ddlam = -3.d0*shearmod - (1.d0 - 2.d0*shearmod*A)*dsigy_ddlam

    dlam = dlam - r / dr_ddlam

    iter = iter + 1

  end do
end subroutine

! ==============================================================================
! SUBROUTINE TRIP_RETURN_MAP
! Return map algorithm for only TRIP

subroutine trip_return_map(sig_eq, sig_eq_tr, shearmod, A, failed)

  use parameters
  implicit none

  ! Output
  real*8, intent(out) :: sig_eq

  ! Input
  real*8, intent(in) :: sig_eq_tr, shearmod, A

  ! Convergence flag
  logical, intent(inout) :: failed

  ! Local
  integer :: iter
  real*8 :: r, dr_dsig

  sig_eq = sig_eq_tr
  iter = 0

  do

    r = sig_eq_tr - (1.d0 + 2.d0*shearmod*A)*sig_eq

    if( abs(r) < param%tolerance * sig_eq_tr ) then
      exit
    end if
    
    if( iter == 100 ) then
      write(7,*) ''
      write(7,*) '*** TRIP RETURN MAP DID NOT CONVERGE'
      failed = .true.
      return
    end if

    dr_dsig = -1.d0 - 2.d0*shearmod*A

    sig_eq = sig_eq - r / dr_dsig
    
    iter = iter + 1

  end do
end subroutine

! ==============================================================================
! SUBROUTINE GET_HARDENING
! Yield strength and hardenings modulus of power law

subroutine get_hardening(sigy, dsigy_ddlam, dlam, sigy0, shearmod, bulkmod, A, lam_tot, failed)

  use parameters
  implicit none

  ! Output
  real*8, intent(out) :: sigy, dsigy_ddlam

  ! Input
  real*8, intent(in) :: dlam, sigy0, shearmod, bulkmod, A, lam_tot

  ! Convergence flag
  logical, intent(inout) :: failed

  ! Local 
  real*8 :: youngs, r, dr_dsigy
  integer :: iter

  youngs = 9.d0*bulkmod*shearmod/(3.d0*bulkmod + shearmod)

  if( A == 0.d0 ) then

    sigy = sigy0*(1.d0 + youngs/sigy0*(lam_tot + dlam))**param%n

  else
    sigy = sigy0
    iter = 0
    do while( .true. )

      r = sigy - sigy0*(1.d0 + youngs/sigy0*(lam_tot + dlam + 2.d0/3.d0*A*sigy))**param%n
    
      if( abs(r) < param%tolerance * 1.d-1 * max(sigy0, 1.d-6) ) then 
        exit
      end if
      
      if( iter == 100 ) then
        write(7,*) ''
        write(7,*) '*** HARDENING DID NOT CONVERGE'
        write(7,*) 'sigy', sigy
        write(7,*) 'sigy0', sigy0
        write(7,*) 'youngs', youngs
        write(7,*) 'lam', lam_tot + dlam
        write(7,*) 'A', A
        failed = .true.
        return
      end if

      dr_dsigy = 1.d0 - 2.d0*A*param%n*youngs/3.d0*(1.d0 + youngs/sigy0*(lam_tot + dlam + 2.d0/3.d0*A*sigy))**(param%n-1.d0)
      sigy = sigy - r / dr_dsigy
      
      iter = iter + 1
    end do
  end if

  dsigy_ddlam = 1.d0 / ( (sigy/sigy0)**(1.d0/param%n-1.d0)/(param%n*youngs) -2.d0/3.d0*A )

end subroutine

! ==============================================================================
! SUBROUTINE READ_PROPS
! Read phase dependent properties from file

subroutine read_props

  use parameters 
  implicit none

  character(len=256):: header
  integer :: npoints, ipoint, iprop
  logical :: exist

  write(6,*) 'READING PROP FILE'
  write(6,*) param%propfile

  inquire(file=trim(param%propfile), exist=exist)
  if( .not. exist ) then 
    write(6,*) 'INPUT FILE NOT FOUND'
    call xit 
  end if

  open(unit=101, file=trim(param%propfile), ACTION='read')
  iprop = 0

  do while( iprop /= 5 )
    read(101,*) header
    read(101,*) npoints

    if( header == "*shear" ) then
      allocate(param%data_shear(npoints,6))
      do ipoint = 1, npoints
         read(101,*) param%data_shear(ipoint,:)
      end do
    elseif( header == "*bulk" ) then
      allocate(param%data_bulk(npoints,6))
      do ipoint = 1, npoints
         read(101,*) param%data_bulk(ipoint,:)
      end do
    elseif( header == "*yield" ) then
      allocate(param%data_yield(npoints,6))
      do ipoint = 1, npoints
         read(101,*) param%data_yield(ipoint,:)
      end do
    elseif( header == "*density" ) then
      allocate(param%data_density(npoints,6))
      do ipoint = 1, npoints
         read(101,*) param%data_density(ipoint,:)
      end do
    elseif( header == "*trip" ) then
      allocate(param%data_trip(npoints,5))
      do ipoint = 1, npoints
         read(101,*) param%data_trip(ipoint,:)
      end do
    else 
      write(7,*) 'INPUTFILE HEADER IS WRONG'
      write(7,*) header
      call xit 
    end if

    iprop = iprop + 1
  end do
  close(unit=101)
end subroutine read_props

! ==============================================================================
! SUBROUTINE INTERP_PROPS
! Calculates phase and temperature dependent properties
! NOTE: Somehow asume-size arrays do not work with Abaqus?????

subroutine interp_prop(prop, temp, table, nphase, ndata)

  implicit none

  ! Input
  real*8, intent(out) :: prop(nphase)
  real*8, intent(in) :: temp
  real*8, dimension(ndata,nphase+1), intent(in) :: table
  integer, intent(in) :: nphase, ndata

  ! Local
  integer :: i, itemp

  itemp = nphase+1

  if( temp <= table(1,itemp) ) then
    prop = table(1,1:nphase)

  elseif( temp >= table(ndata,itemp) ) then
    prop = table(ndata,1:nphase)
    
  else
    do i = 1, ndata - 1
      if ( temp >= table(i,itemp) .and. temp < table(i+1,itemp)) then
        prop = table(i,1:nphase) &
          + (temp - table(i,itemp))/(table(i+1,itemp) - table(i,itemp)) &
          * (table(i+1,1:nphase) - table(i,1:nphase))
        exit
      end if
    end do
  end if

end subroutine interp_prop

! ==============================================================================
! SUBROUTINE INTERP_PROPS
! Calculates phase and temperature dependent properties

! subroutine interp_props_old(values, x)

!   use parameters
!   implicit none

!   ! Input
!   real*8, dimension(24), intent(out) :: values
!   real*8, intent(in) :: x

!   ! Local
!   integer :: ndata, i

!   ndata = size(param%data_mech, 1)

!   if( x <= param%data_mech(1,1) ) then
!     values = param%data_mech(1,2:25)

!   elseif( x >= param%data_mech(ndata,1) ) then
!     values = param%data_mech(ndata,2:25)

!   else
!     do i = 1, ndata - 1
!       if ( x >= param%data_mech(i,1) .and. x < param%data_mech(i+1,1)) then
!         values = param%data_mech(i,2:25) &
!           + (x - param%data_mech(i,1))/(param%data_mech(i+1,1) - param%data_mech(i,1)) &
!           * (param%data_mech(i+1,2:25) - param%data_mech(i,2:25))
!         exit
!       end if
!     end do
!   end if

! end subroutine

