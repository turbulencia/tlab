#include "dns_error.h"
!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/05 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   Compute gammas - needed for conditional/intrinsic 
!#   averaging of 1d-vertical-statistics [cf. Pope. p.170].
!#   
!#   The gammas describe the share of different partitions in the 
!#   whole simulation domain:
!#
!#      gamma_1: partition solid (with interfaces) 
!#               [solid and interface partitions are described 
!#               with ones in the eps - indicator field]
!#      gamma_0: partition fluid (without interfaces) 
!#               [fluid partitions are described with zeros in 
!#               in the eps - indicator field]
!#
!#   And the physically correct partitions (based on volume fractions)
!#   [interface points are splitted in half for solid and fluid]:
!#      {the following is not valid at edges and corners!
!#       --> special treatment of these points (e.g. for gamma_s):
!#           - plain interface : 1/2
!#           - edges           : 1/4
!#           - corners         : 1/8}
!#      gamma_f: partition fluid
!#      gamma_s: partition solid
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

subroutine IBM_AVG_GAMMA(gamma_0, gamma_1, gamma_f, gamma_s, eps, tmp1, tmp2)

  use TLAB_VARS,      only : imax, jmax, kmax, g, area
  use TLAB_CONSTANTS, only : efile, wi, wp
  use TLAB_PROCS,     only : TLAB_STOP, TLAB_WRITE_ASCII
  use AVGS,           only : AVG_IK_V
  use IBM_VARS,       only : dy, facu, facl
#ifdef USE_MPI
  use MPI
  use TLAB_MPI_VARS,  only : ims_err 
#endif

  implicit none

  real(wp), dimension(     jmax     ), intent(out  ) :: gamma_0, gamma_1
  real(wp), dimension(     jmax     ), intent(out  ) :: gamma_f, gamma_s
  real(wp), dimension(imax,jmax,kmax), intent(in   ) :: eps
  real(wp), dimension(imax,jmax,kmax), intent(inout) :: tmp1
  real(wp), dimension(     jmax,   2), intent(inout) :: tmp2

  target                                             :: tmp2
  real(wp), dimension(:),              pointer       :: wrk1d, res

  integer(wi)                                        :: i,j,k
  real(wp)                                           :: dummy, check
  
  ! ================================================================== !
  ! solid + interface points are filled with zeros
  tmp1(:,:,:) = 0.0_wp; tmp2(:,:) = 0.0_wp
  wrk1d => tmp2(:,1); res => tmp2(:,2)

  ! horizontal average - compute gamma_1
  CALL AVG_IK_V(imax,jmax,kmax, jmax, eps,  g(1)%jac, g(3)%jac, gamma_1, wrk1d, area)
  
  ! horizontal average - compute gamma_0
  tmp1(:,:,:) = (1.0_wp - eps(:,:,:))
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac, g(3)%jac, gamma_0, wrk1d, area)

  ! take boundaries, edges and corner into account
  tmp1(:,:,:) = eps(:,:,:)

  ! take grid stretching in vertical direction into account
  do j = 1, jmax-1   ! dy[jmax-1] - vertical spacing
    dy(j) = g(2)%nodes(j+1) - g(2)%nodes(j) 
  end do
  do j = 1, jmax-2   ! fac[jmax-2]
    facl(j) = dy(j  ) / ( dy(j) + dy(j+1) ) ! lower obj
    facu(j) = dy(j+1) / ( dy(j) + dy(j+1) ) ! upper obj
  end do

  ! in x-direction - inner points - boundary points
  do k = 1, kmax   
    do j = 1, jmax 
      i = 1
      if (        eps(i,j,k) == 0 .and. eps(i+1,j,k) == 1 ) then
        tmp1(i+1,j,k) = 0.5_wp * tmp1(i+1,j,k)
      end if
      do i = 2, imax-1
        if (      eps(i,j,k) == 0 .and. eps(i+1,j,k) == 1 ) then
          tmp1(i+1,j,k) = 0.5_wp * tmp1(i+1,j,k)
        else if ( eps(i,j,k) == 0 .and. eps(i-1,j,k) == 1 ) then
          tmp1(i-1,j,k) = 0.5_wp * tmp1(i-1,j,k)
        end if
      end do 
      i = imax
      if (        eps(i,j,k) == 0 .and. eps(i-1,j,k) == 1 ) then
        tmp1(i-1,j,k) = 0.5_wp * tmp1(i-1,j,k)
      end if
    end do

    ! in y-direction - inner points - boundary points 
    ! (without considering stretching on the ground)
    do i = 1, imax   
      j = 1        
      if (        eps(i,j,k) == 0 .and. eps(i,j+1,k) == 1 ) then
        tmp1(i,j+1,k) = 0.5_wp * tmp1(i,j+1,k)   
      end if 
      do j = 2, jmax-1
        if (      eps(i,j,k) == 0 .and. eps(i,j+1,k) == 1 ) then
          if ( j == 2 .or. j == jmax-1 ) then
            tmp1(i,j+1,k) = 0.5_wp * tmp1(i,j+1,k) 
          else
            tmp1(i,j+1,k) = facu(j) * tmp1(i,j+1,k)
          end if
          ! tmp1(i,j+1,k) = 0.5_wp * tmp1(i,j+1,k)
        else if ( eps(i,j,k) == 0 .and. eps(i,j-1,k) == 1 ) then
          if ( j == 2 .or. j == jmax-1 ) then
            tmp1(i,j-1,k) = 0.5_wp * tmp1(i,j-1,k)
          else
            tmp1(i,j-1,k) = facl(j-2) * tmp1(i,j-1,k)
          end if
          ! tmp1(i,j-1,k) = 0.5_wp * tmp1(i,j-1,k) 
        end if
      end do 
      j = jmax
      if (        eps(i,j,k) == 0 .and. eps(i,j-1,k) == 1 ) then
        tmp1(i,j-1,k) = 0.5_wp * tmp1(i,j-1,k)
      end if
    end do  
  end do

  ! in z-direction - inner points - boundary points
  do i = 1, imax   
    do j = 1, jmax 
      k = 1        
      if (        eps(i,j,k) == 0 .and. eps(i,j,k+1) == 1 ) then
        tmp1(i,j,k+1) = 0.5_wp * tmp1(i,j,k+1)
      end if 
      do k = 2, kmax-1
        if (      eps(i,j,k) == 0 .and. eps(i,j,k+1) == 1 ) then
          tmp1(i,j,k+1) = 0.5_wp * tmp1(i,j,k+1)
        else if ( eps(i,j,k) == 0 .and. eps(i,j,k-1) == 1 ) then
          tmp1(i,j,k-1) = 0.5_wp * tmp1(i,j,k-1)
        end if
      end do 
      k = kmax
      if (        eps(i,j,k) == 0 .and. eps(i,j,k-1) == 1 ) then
        tmp1(i,j,k-1) = 0.5_wp * tmp1(i,j,k-1)
      end if
    end do
  end do

  ! horizontal average - compute gamma_s
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac, g(3)%jac, gamma_s, wrk1d, area)

  ! horizontal average - compute gamma_f
  tmp1(:,:,:) = (1.0_wp - tmp1(:,:,:))
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac, g(3)%jac, gamma_f, wrk1d, area)

  ! constistency check of partitioning
  do j = 1, jmax 
    if ( gamma_f(j) /= 1.0_wp ) then
      res(j) = (gamma_f(j) - gamma_0(j)) + (gamma_s(j) - gamma_1(j)) 
    else
      res(j) = 0.0_wp
    end if
  end do
  check = sum(res)
#ifdef USE_MPI
  dummy = check
  CALL MPI_ALLREDUCE(dummy, check, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
  if ( check > 1.0e-10_wp ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_AVG_GAMMA. Sum of gammas is not correct.')
    call TLAB_STOP(DNS_ERROR_IBM_GAMMA)
  end if

  nullify(wrk1d, res)

  return
end subroutine IBM_AVG_GAMMA

!########################################################################

subroutine IBM_AVG_SCAL_BCS(is, scalv_bcs)
  
  use IBM_VARS,       only : ibm_objup, max_height_objlo, max_height_objup
  use IBM_VARS,       only : ibmscaljmin, ibmscaljmax
  use TLAB_VARS,      only : jmax
  use TLAB_CONSTANTS, only : wi, wp

  implicit none

  integer(wi),               intent(in ) :: is
  real(wp), dimension(jmax), intent(out) :: scalv_bcs
  
  integer(wi)                            :: j

  ! ================================================================== !
  ! write out scalar boundary values applied in solids (vertical)
  ! assuming homogenous temperature in solids, otherwise change here
  
  scalv_bcs = 0.0_wp

  do j = 1, int(max_height_objlo, wi)
    scalv_bcs(j) = ibmscaljmin(is)
  end do
  if ( ibm_objup ) then
    do j = jmax-int(max_height_objup, wi),jmax
      scalv_bcs(j) = ibmscaljmax(is)
    end do
  end if

  return
end subroutine IBM_AVG_SCAL_BCS