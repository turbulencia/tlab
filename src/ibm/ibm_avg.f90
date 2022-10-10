#include "types.h"
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
  use TLAB_CONSTANTS, only : efile
  use TLAB_PROCS,     only : TLAB_STOP, TLAB_WRITE_ASCII
#ifdef USE_MPI
  use MPI
  use TLAB_MPI_VARS,  only : ims_err 
#endif
  use IO_FIELDS

  implicit none

  TREAL, dimension(     jmax     ), intent(out  ) :: gamma_0, gamma_1
  TREAL, dimension(     jmax     ), intent(out  ) :: gamma_f, gamma_s
  TREAL, dimension(imax,jmax,kmax), intent(in   ) :: eps
  TREAL, dimension(imax,jmax,kmax), intent(inout) :: tmp1
  TREAL, dimension(jmax,2),         intent(inout) :: tmp2

  target                                          :: tmp2
  TREAL, dimension(:),              pointer       :: wrk1d, res
  
  TINTEGER                                        :: i,j,k
  TREAL                                           :: dummy, check
  
  TREAL, dimension(:,:,:),allocatable :: tmp
  allocate(tmp(imax,jmax,kmax))

  ! ================================================================== !
  ! solid + interface points are filled with zeros
  tmp1(:,:,:) = C_0_R; tmp2(:,:) = C_0_R
  wrk1d => tmp2(:,1); res => tmp2(:,2)

  ! horizontal average - compute gamma_1
  CALL AVG_IK_V(imax,jmax,kmax, jmax, eps,  g(1)%jac, g(3)%jac, gamma_1, wrk1d, area)
  
  ! horizontal average - compute gamma_0
  tmp1(:,:,:) = (C_1_R - eps(:,:,:))
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac, g(3)%jac, gamma_0, wrk1d, area)

  ! take boundaries, edges and corner into account
  tmp1(:,:,:) = eps(:,:,:)

  ! in x-direction - inner points - boundary points
  do k = 1, kmax   
    do j = 1, jmax 
      i = 1
      if (        eps(i,j,k) == 0 .and. eps(i+1,j,k) == 1 ) then
        tmp1(i+1,j,k) = C_05_R * tmp1(i+1,j,k)
      end if
      do i = 2, imax-1
        if (      eps(i,j,k) == 0 .and. eps(i+1,j,k) == 1 ) then
          tmp1(i+1,j,k) = C_05_R * tmp1(i+1,j,k)
        else if ( eps(i,j,k) == 0 .and. eps(i-1,j,k) == 1 ) then
          tmp1(i-1,j,k) = C_05_R * tmp1(i-1,j,k)
        end if
      end do 
      i = imax
      if (        eps(i,j,k) == 0 .and. eps(i-1,j,k) == 1 ) then
        tmp1(i-1,j,k) = C_05_R * tmp1(i-1,j,k)
      end if
    end do

    ! in y-direction - inner points - boundary points
    do i = 1, imax   
      j = 1        
      if (        eps(i,j,k) == 0 .and. eps(i,j+1,k) == 1 ) then
        tmp1(i,j+1,k) = C_05_R * tmp1(i,j+1,k)
      end if 
      do j = 2, jmax-1
        if (      eps(i,j,k) == 0 .and. eps(i,j+1,k) == 1 ) then
          tmp1(i,j+1,k) = C_05_R * tmp1(i,j+1,k)
        else if ( eps(i,j,k) == 0 .and. eps(i,j-1,k) == 1 ) then
          tmp1(i,j-1,k) = C_05_R * tmp1(i,j-1,k)
        end if
      end do 
      j = jmax
      if (        eps(i,j,k) == 0 .and. eps(i,j-1,k) == 1 ) then
        tmp1(i,j-1,k) = C_05_R * tmp1(i,j-1,k)
      end if
    end do  
  end do

  ! in z-direction - inner points - boundary points
  do i = 1, imax   
    do j = 1, jmax 
      k = 1        
      if (        eps(i,j,k) == 0 .and. eps(i,j,k+1) == 1 ) then
        tmp1(i,j,k+1) = C_05_R * tmp1(i,j,k+1)
      end if 
      do k = 2, kmax-1
        if (      eps(i,j,k) == 0 .and. eps(i,j,k+1) == 1 ) then
          tmp1(i,j,k+1) = C_05_R * tmp1(i,j,k+1)
        else if ( eps(i,j,k) == 0 .and. eps(i,j,k-1) == 1 ) then
          tmp1(i,j,k-1) = C_05_R * tmp1(i,j,k-1)
        end if
      end do 
      k = kmax
      if (        eps(i,j,k) == 0 .and. eps(i,j,k-1) == 1 ) then
        tmp1(i,j,k-1) = C_05_R * tmp1(i,j,k-1)
      end if
    end do
  end do
  
  CALL IO_WRITE_FIELDS('flow.10', IO_FLOW, imax,jmax,kmax, 1, tmp1, tmp)

  ! horizontal average - compute gamma_s
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac, g(3)%jac, gamma_s, wrk1d, area)

  ! horizontal average - compute gamma_f
  tmp1(:,:,:) = (C_1_R - tmp1(:,:,:))
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac, g(3)%jac, gamma_f, wrk1d, area)

  ! constistency check of partitioning
  do j = 1, jmax 
    if ( gamma_f(j) /= C_1_R ) then
      res(j) = (gamma_f(j) - gamma_0(j)) + (gamma_s(j) - gamma_1(j)) 
    else
      res(j) = C_0_R
    end if
  end do
  check = sum(res)
#ifdef USE_MPI
  dummy = check
  CALL MPI_ALLREDUCE(dummy, check, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
  if ( check > C_1EM10_R ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_AVG_GAMMA. Sum of gammas is not correct.')
    call TLAB_STOP(DNS_ERROR_IBM_GAMMA)
  end if

  nullify(wrk1d, res)

  return
end subroutine IBM_AVG_GAMMA

!########################################################################

subroutine IBM_AVG_SCAL_BCS(is, scalv_bcs)
  
  use IBM_VARS,  only : ibm_objup, max_height_objlo, max_height_objup
  use IBM_VARS,  only : ibmscaljmin, ibmscaljmax
  use TLAB_VARS, only : jmax

  implicit none

  TINTEGER,               intent(in ) :: is
  TREAL, dimension(jmax), intent(out) :: scalv_bcs
  
  TINTEGER                            :: j

  ! ================================================================== !
  ! write out scalar boundary values applied in solids (vertical)
  ! assuming homogenous temperature in solids, otherwise change here
  
  scalv_bcs = C_0_R

  do j = 1, int(max_height_objlo)
    scalv_bcs(j) = ibmscaljmin(is)
  end do
  if ( ibm_objup ) then
    do j = jmax-int(max_height_objup),jmax
      scalv_bcs(j) = ibmscaljmax(is)
    end do
  end if

  return
end subroutine IBM_AVG_SCAL_BCS