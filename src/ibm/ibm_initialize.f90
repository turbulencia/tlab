#include "dns_error.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#
!#   called once in dns_main.f90 before time integration starts,
!#   all relevant geometry informations needed for the IBM are 
!#   available after calling this routine
!#    
!#   options for geometry generation:
!#      0. use existing eps field (ibm_restart==.true.)
!#      1. generate eps field outside tlab and use as existing eps field with
!#         (ibm_restart==.true.)
!#      2. write own routine to generate geometry 
!#         (cf. IBM_GENERATE_GEOMETRY_XBARS, ibm_restart==.false.)
!# 
!########################################################################
!# ARGUMENTS 
!#
!# 
!########################################################################
!# REQUIREMENTS
!#
!# 
!########################################################################

subroutine IBM_INITIALIZE_GEOMETRY(txc, wrk3d)  
  
  use IBM_VARS
  use TLAB_VARS,      only : isize_field, inb_txc
  use TLAB_VARS,      only : stagger_on
  use TLAB_CONSTANTS, only : efile, wp
  use IO_FIELDS
  use TLAB_PROCS

  implicit none

  real(wp), dimension(isize_field,inb_txc), intent(inout) :: txc
  real(wp), dimension(isize_field),         intent(inout) :: wrk3d

  target                                                  :: txc

  real(wp), dimension(:), pointer                         :: epsi, epsj, epsk
  real(wp), dimension(:), pointer                         :: tmp1, tmp2
#ifdef IBM_DEBUG
  real(wp), dimension(:), pointer                         :: tmp3
#endif
  logical :: flag_epsp
  ! ================================================================== !
  ! assigning pointer to scratch
  txc = 0.0_wp; epsi => txc(:,1); epsj => txc(:,2); epsk => txc(:,3)
  tmp1 => txc(:,4); tmp2 => txc(:,5)

  ! eps field (read/create)
  if ( ibm_restart ) then
    flag_epsp = .false.
    call IBM_IO_READ(wrk3d, flag_epsp)
  else if ( xbars_geo%name == 'xbars' ) then
    call IBM_GENERATE_GEOMETRY_XBARS(wrk3d)
  else
    call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY no objects in flow.')
    call TLAB_STOP(DNS_ERROR_IBM_MISS_GEO)
  end if

  ! transpose eps (epsi, epsj, epsk)
  call IBM_GEOMETRY_TRANSPOSE(epsi, epsj, epsk, wrk3d)

  ! generate relevant geometry fields for IBM routines (nobi, nobj, nobk)
  call IBM_GENERATE_GEOMETRY(epsi, epsj, epsk)

  ! verify geometry
  call IBM_VERIFY_GEOMETRY()

  ! generate epsp on pressure mesh
  if ( stagger_on ) then
    if ( ibm_restart ) then
      flag_epsp = .true.
      call IBM_IO_READ(wrk3d, flag_epsp)
    else if (.not. (xbars_geo%name == 'xbars') ) then
      call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY epsp field is missing.')
      call TLAB_STOP(DNS_ERROR_IBM_MISS_GEO)
    end if
  end if

  ! compute gammas for conditional averages 
  CALL IBM_AVG_GAMMA(gamma_0, gamma_1, eps, tmp1)
  tmp1(:) = 0.0_wp

  ! check idle procs
#ifdef USE_MPI
  call IBM_CHECK_PROCS(epsi, epsj, epsk)
#else   
  ! in case of serial mode: one task with full domain, no idle procs
  ims_pro_ibm_x = .true.; ims_pro_ibm_y = .true.; ims_pro_ibm_z = .true.
#endif

#ifdef IBM_DEBUG
  ! io of all geometry fields in debugging mode 
  tmp3 => txc(:,6); tmp3(:) = 0.0_wp
  call IBM_GEOMETRY_DEBUG_IO(epsi, epsj, epsk, tmp1, tmp2, tmp3)
  nullify(tmp3)
#endif

  ! switch to true in routines where the IBM is needed
  ibm_burgers = .false.; ibm_partial = .false.
  
  ! disassociate pointers
  nullify(epsi, epsj, epsk)
  nullify(tmp1, tmp2)
  wrk3d(:) = 0.0_wp

  if (stagger_on)  then
    write(*,*)'sum,min,max of eps  :', sum(eps), minval(eps), maxval(eps)
    write(*,*)'sum,min,max of epsp :', sum(epsp), minval(epsp), maxval(epsp)
    call IBM_STAGGER_GEOMETRY(eps, wrk3d)
    write(*,*) '##########################################'
    write(*,*) 'difference sum(eps  -epsp   ):', sum(eps-epsp)
    write(*,*) 'difference sum(eps  -routine):', sum(eps-wrk3d)
    write(*,*) 'difference sum(epsp -routine):', sum(epsp-wrk3d)
    wrk3d(:) = 0.0_wp
  end if

  return
end subroutine IBM_INITIALIZE_GEOMETRY

!########################################################################

subroutine IBM_IO_READ(wrk3d, flag_epsp)

  use IBM_VARS
  use TLAB_VARS,      only : imax,jmax,kmax, isize_field
  use TLAB_CONSTANTS, only : wp, wi
  use IO_FIELDS
  
  implicit none
  
  real(wp), dimension(isize_field), intent(inout) :: wrk3d
  logical,                          intent(inout) :: flag_epsp
  
  character(len=32)                               :: name
  
    ! ================================================================== !
  wrk3d(:) = 0.0_wp
  select case( ibm_io )
    case ( IBM_IO_REAL )
      if (flag_epsp) then
        name = epsp_name_real
        call IO_READ_FIELDS(name, IO_FLOW, imax,jmax,kmax, 1, 0, epsp)
      else
        name = eps_name_real
        call IO_READ_FIELDS(name, IO_FLOW, imax,jmax,kmax, 1, 0, eps)
      end if
    case ( IBM_IO_INT  )
      call IBM_IO_READ_INT_GEOMETRY(wrk3d, flag_epsp)
    case ( IBM_IO_BIT  )
      call IBM_IO_READ_BIT_GEOMETRY(wrk3d, flag_epsp)
  end select 
  
  return
end subroutine IBM_IO_READ

!########################################################################

subroutine IBM_IO_WRITE(wrk3d, flag_epsp)

  use IBM_VARS
  use TLAB_VARS,      only : imax,jmax,kmax, isize_field
  use TLAB_CONSTANTS, only : wp, wi
  use IO_FIELDS
  
  implicit none
  
  real(wp), dimension(isize_field), intent(inout) :: wrk3d
  logical,                          intent(inout) :: flag_epsp
  
  character(len=32)                               :: name
  
    ! ================================================================== !
  wrk3d(:) = 0.0_wp
  select case( ibm_io )
    case ( IBM_IO_REAL )
      if (flag_epsp) then
        name = epsp_name_real
        call IO_WRITE_FIELDS(name, IO_FLOW, imax,jmax,kmax, 1, epsp)
      else
        name = eps_name_real
        call IO_WRITE_FIELDS(name, IO_FLOW, imax,jmax,kmax, 1, eps)
      end if
    case ( IBM_IO_INT  )
      call IBM_IO_WRITE_INT_GEOMETRY(wrk3d, flag_epsp)
    case ( IBM_IO_BIT  )
      call IBM_IO_WRITE_BIT_GEOMETRY(wrk3d, flag_epsp)
  end select 
  
  return
end subroutine IBM_IO_WRITE

!########################################################################

subroutine IBM_STAGGER_GEOMETRY(eps, epsp)

  use TLAB_VARS,      only : imax, jmax, kmax
  use TLAB_CONSTANTS, only : wp, wi

  implicit none

  real(wp), dimension(imax,jmax,kmax), intent(in ) :: eps
  real(wp), dimension(imax,jmax,kmax), intent(out) :: epsp

  integer(wi)                                      :: i,j,k

  ! ================================================================== !
   
  ! horizontal staggering - move right interfaces one step to the left,
  ! don't touch vertical direction since no staggering is applied here
 
  epsp(:,:,:) = eps(:,:,:)
  
  do k = 1, kmax   
    do j = 1, jmax 
      do i = 2, imax
        if ( eps(i,j,k) == 0 .and. eps(i-1,j,k) == 1 ) then
          epsp(i-1,j,k) = 0.0_wp 
        end if
      end do 
    end do
  end do
  
  do i = 1, imax   
    do j = 1, jmax 
      do k = 2, kmax
        if ( eps(i,j,k) == 0 .and. eps(i,j,k-1) == 1 ) then
          epsp(i,j,k-1) = 0.0_wp 
        end if
      end do 
    end do
  end do

  return
end subroutine IBM_STAGGER_GEOMETRY

!########################################################################