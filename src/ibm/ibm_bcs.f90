!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   apply IBM boundary conditions to flow/scalar fields
!#  
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

subroutine IBM_INITIALIZE_SCAL(isbcs, s)
  
  use IBM_VARS,       only : eps, ibmscaljmin, ibmscaljmax, scal_bcs
  use TLAB_VARS,      only : imax, jmax, isize_field, inb_scal
  use TLAB_CONSTANTS, only : wi, wp

  implicit none

  integer(wi),                                  intent(in   ) :: isbcs
  real(wp),    dimension(isize_field,inb_scal), intent(inout) :: s

  integer(wi)                                                 :: is

! ================================================================== !
! get scalar dirichlet boundary values of ini scalar field
! (assuming always homogenous horizontal temperature on upper/lower boundaries)
  do is = 1, inb_scal
    ibmscaljmin(is) = s(1,                is)
    ibmscaljmax(is) = s(imax*(jmax-1) + 1,is)
  end do

! set scalar values in solid to zero
  call IBM_BCS_FIELD_COMBINED(1,s)

! apply ibmscaljmin, ibmscaljmax on scalar field(s)
  do is = 1, inb_scal
    call IBM_BCS_SCAL(is,s(:,is),eps(:))
  end do

! write out scalar boundary values applied in solids
  if ( isbcs > 0 ) then
    do is = 1, inb_scal
      call IBM_AVG_SCAL_BCS(is, scal_bcs(:,is))
    end do
  end if

  return
end subroutine IBM_INITIALIZE_SCAL

!########################################################################

subroutine IBM_BCS_SCAL(is,s,eps)
  
  use IBM_VARS,       only : ibmscaljmin, ibmscaljmax, ibm_objup, max_height_objup
  use TLAB_VARS,      only : imax, jmax, kmax
  use TLAB_CONSTANTS, only : wi, wp

  implicit none

  integer(wi),                            intent(in   ) :: is
  real(wp),    dimension(imax,jmax,kmax), intent(inout) :: s
  real(wp),    dimension(imax,jmax,kmax), intent(in   ) :: eps

  integer(wi)                                           :: j

! ================================================================== !
! default, set only scalar value on lower boundary
! if objects also on upper boundary present, 
! values needs to be overwritten, cf. next loop
  s(:,:,:) = s(:,:,:) + eps(:,:,:) * ibmscaljmin(is) 

! in case of objects on upper boundary, set different temperature here
  if ( ibm_objup ) then
    do j = jmax-int(max_height_objup),jmax
      s(:,j,:) = (1.0_wp - eps(:,j,:)) *  s(:,j,:) + eps(:,j,:) * ibmscaljmax(is) 
    end do
  end if 

  return
end subroutine IBM_BCS_SCAL

!########################################################################

subroutine IBM_BCS_FIELD_COMBINED(is,fld)
  
  use TLAB_VARS,      only : isize_field, inb_flow, inb_scal
  use TLAB_CONSTANTS, only : wi, wp

  implicit none
  
  integer(wi),                        intent(in   ) :: is
  real(wp), dimension(isize_field,*), intent(inout) :: fld

  integer(wi)                                       :: i

  ! ================================================================== !
  ! apply IBM BCs on many fields
  if ( is == 0 ) then
    do i = 1, inb_flow
      call IBM_BCS_FIELD(fld(1,i))
    end do
  elseif ( is == 1 ) then
    do i = 1, inb_scal
      call IBM_BCS_FIELD(fld(1,i))
    end do
  end if

  return
end subroutine IBM_BCS_FIELD_COMBINED

!########################################################################

subroutine IBM_BCS_FIELD(fld)
  
  use IBM_VARS,       only : eps
  use TLAB_VARS,      only : isize_field
  use TLAB_CONSTANTS, only : wi, wp

  implicit none
  
  real(wp), dimension(isize_field), intent(inout) :: fld

  integer(wi)                                     :: i

  ! ================================================================== !
  ! apply IBM BCs on scalar/flow fields
  do i = 1, isize_field
    fld(i) = (1.0_wp - eps(i)) * fld(i)  
  end do

  return
end subroutine IBM_BCS_FIELD

!########################################################################

subroutine IBM_BCS_FIELD_STAGGER(fld)
  
  use IBM_VARS,       only : epsp
  use TLAB_VARS,      only : isize_field
  use TLAB_CONSTANTS, only : wi, wp

  implicit none
  
  real(wp), dimension(isize_field), intent(inout) :: fld

  integer(wi)                                     :: i

  ! ================================================================== !
  ! apply IBM BCs on scalar/flow fields
  do i = 1, isize_field
    fld(i) = (1.0_wp - epsp(i)) * fld(i)  
  end do

  return
end subroutine IBM_BCS_FIELD_STAGGER

!########################################################################

subroutine IBM_BCS_FIELD_INV(fld,tmp) ! not used so far
  
  use IBM_VARS,       only : eps
  use TLAB_VARS,      only : isize_field
  use TLAB_CONSTANTS, only : wi, wp

  implicit none
  
  real(wp), dimension(isize_field), intent(in)  :: fld
  real(wp), dimension(isize_field), intent(out) :: tmp

  integer(wi)                                   :: i

  ! ================================================================== !
  ! apply inverse IBM BCs on fields -- only BCs in solid left, fluid regions are zero
  do i = 1, isize_field
    tmp(i) = eps(i) * fld(i)  
  end do

  return
end subroutine IBM_BCS_FIELD_INV