#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2021/XX/XX - J.K. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!#
!#
!#
!#
!########################################################################
subroutine FI_CHANNEL_CFR_FORCING(u, v, h1, h2, tmp1, wrk1d, wrk2d, wrk3d)

  use TLAB_VARS, only : g, imax,jmax,kmax, isize_field,isize_wrk2d
  use TLAB_VARS, only : ubulk, ubulk_parabolic

  implicit none

#include "integers.h"

  TREAL, dimension(isize_field), intent(inout) :: u     ! streamwise velocity
  TREAL, dimension(isize_field), intent(inout) :: v     ! streamwise velocity
  TREAL, dimension(isize_field), intent(inout) :: h1    
  TREAL, dimension(isize_field), intent(inout) :: h2    
  TREAL, dimension(isize_field), intent(inout) :: tmp1
  TREAL, dimension(jmax),        intent(inout) :: wrk1d
  TREAL, dimension(isize_wrk2d), intent(inout) :: wrk2d
  TREAL, dimension(isize_field), intent(inout) :: wrk3d

! -----------------------------------------------------------------------

  TINTEGER                                     :: ij
  TINTEGER, dimension(2,2)                     :: bcs
  TREAL                                        :: delta_ub, wrot

! #######################################################################

  call FI_CHANNEL_UBULK(u, wrk1d, wrk3d) 

  delta_ub = - (ubulk_parabolic - ubulk)

  ! write(*,*)  ubulk_parabolic, ubulk, delta_ub ! debug

  do ij = 1, isize_field
    h1(ij) = h1(ij) - delta_ub
  end do
  
  return
end subroutine FI_CHANNEL_CFR_FORCING
!########################################################################
subroutine FI_CHANNEL_CPG_FORCING(u, v, h1, h2, tmp1, wrk1d, wrk2d, wrk3d)

  use TLAB_VARS, only : g, imax,jmax,kmax, isize_field,isize_wrk2d
  use TLAB_VARS, only : ubulk, ubulk_parabolic, visc, f_cpg

  implicit none

#include "integers.h"

  TREAL, dimension(isize_field), intent(inout) :: u     ! streamwise velocity
  TREAL, dimension(isize_field), intent(inout) :: v     ! streamwise velocity
  TREAL, dimension(isize_field), intent(inout) :: h1    
  TREAL, dimension(isize_field), intent(inout) :: h2    
  TREAL, dimension(isize_field), intent(inout) :: tmp1
  TREAL, dimension(jmax),        intent(inout) :: wrk1d
  TREAL, dimension(isize_wrk2d), intent(inout) :: wrk2d
  TREAL, dimension(isize_field), intent(inout) :: wrk3d

! -----------------------------------------------------------------------

  TINTEGER                                     :: ij
  ! TINTEGER, dimension(2,2)                     :: bcs
  TREAL                                        :: delta_ub, wrot

! #######################################################################

  call FI_CHANNEL_UBULK(u, wrk1d, wrk3d) 

  ! write(*,*) reynolds, visc, f_cpg

  do ij = 1, isize_field
    h1(ij) = h1(ij) + f_cpg
  end do

  ! wrot = 0.12

  ! bcs = 0
  ! ! call OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, tmp1,  wrk3d, wrk2d,wrk3d)
  ! do ij = 1, isize_field
  !   h1(ij) = h1(ij) - wrot * v(ij)
  ! end do

  ! ! CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), u, tmp1,  wrk3d, wrk2d,wrk3d)
  ! do ij = 1, isize_field
  !   h2(ij) = h2(ij) + wrot * u(ij)
  ! end do


  return
end subroutine FI_CHANNEL_CPG_FORCING
!########################################################################
subroutine FI_CHANNEL_UBULK(u,wrk1d,wrk3d)

  use TLAB_VARS, only : imax,jmax,kmax, isize_field
  use TLAB_VARS, only : g, area
  use TLAB_VARS, only : ubulk

  implicit none

#include "integers.h"

  TREAL, dimension(isize_field), intent(inout) :: u
  TREAL, dimension(jmax),        intent(inout) :: wrk1d
  TREAL, dimension(isize_field), intent(inout) :: wrk3d

! -----------------------------------------------------------------------
  TREAL, dimension(jmax)                       :: uxz
  TREAL                                        :: SIMPSON_NU

! #######################################################################

  call AVG_IK_V(imax,jmax,kmax, jmax, u, g(1)%jac,g(3)%jac, uxz, wrk1d, area)

  ubulk = SIMPSON_NU(jmax, uxz, g(2)%nodes)

  return
end subroutine FI_CHANNEL_UBULK

!########################################################################

subroutine FI_CHANNEL_UBULK_INITIALIZE()

  use TLAB_VARS,      only : g, qbg, ubulk_parabolic
  use TLAB_CONSTANTS, only : efile
  use TLAB_PROCS

  implicit none

! -----------------------------------------------------------------------

  ! laminar streamwise bulk velocity, make sure that the initial parabolic 
  ! velocity profile is correct (u(y=0)=u(y=Ly)=0)

  if(qbg(1)%type == PROFILE_PARABOLIC) then
    ubulk_parabolic = (C_2_R/C_3_R) * g(2)%nodes(g(2)%size) * qbg(1)%delta 
  else 
    call TLAB_WRITE_ASCII(efile, 'Analytical ubulk cannot be computed. Check initial velocity profile.')
    call TLAB_STOP(DNS_ERROR_UNDEVELOP)
  end if

  return
end subroutine FI_CHANNEL_UBULK_INITIALIZE

!########################################################################