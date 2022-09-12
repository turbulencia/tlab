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
!#      gamma_0: partition solid (without interfaces) 
!#               [solid and interface partitions are described 
!#               with zeros in the (1-eps) - field]
!#      gamma_i: partition interfaces
!#      gamma_1: partition fluid (without interfaces) 
!#               [fluid partitions are described with ones in 
!#               in the (1-eps) - field]
!#
!#   And the physically correct partitions
!#   [interface points are splitted in half to solid and fluid]:
!#      {the following is not valid at edges and corners!
!#       --> special treatment of these points}
!#      gamma_f: partition fluid (N_F = N_1 + 1/2*N_i) 
!#      gamma_s: partition solid (N_S = N_0 + 1/2*N_i)
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

subroutine IBM_AVG_GAMMA(gamma_0, gamma_i, gamma_1, gamma_f, gamma_s, eps, tmp1, tmp2, tmp3)

  use TLAB_VARS,      only : imax, jmax, kmax, g, area
  use TLAB_CONSTANTS, only : efile
  use TLAB_PROCS,     only : TLAB_STOP, TLAB_WRITE_ASCII
#ifdef USE_MPI
  use MPI
  use TLAB_MPI_VARS,  only : ims_err 
#endif

  implicit none

  TREAL, dimension(     jmax     ), intent(out  ) :: gamma_0, gamma_i, gamma_1
  TREAL, dimension(     jmax     ), intent(out  ) :: gamma_f, gamma_s
  TREAL, dimension(imax,jmax,kmax), intent(in   ) :: eps
  TREAL, dimension(imax,jmax,kmax), intent(inout) :: tmp1, tmp2
  TREAL, dimension(jmax,3),         intent(inout) :: tmp3

  target                                          :: tmp3
  TREAL, dimension(:),              pointer       :: wrk1d, res_1, res_2
  
  TINTEGER                                        :: i,j,k
  TREAL                                           :: dummy, check_1, check_2

  ! ================================================================== !
  ! solid + interface points are filled with zeros
  tmp1(:,:,:) = C_0_R; tmp2(:,:,:) = C_0_R; tmp3(:,:) = C_0_R
  wrk1d => tmp3(:,1); res_1 => tmp3(:,2); res_2 => tmp3(:,3)

  ! horizontal average - compute gamma_1
  tmp1(:,:,:) = (C_1_R - eps(:,:,:))
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac, g(3)%jac, gamma_1, tmp3, area)
  tmp1(:,:,:) = eps(:,:,:); tmp2(:,:,:) = eps(:,:,:)

  ! in x-direction - inner points - boundary points
  do k = 1, kmax   
    do j = 1, jmax 
      i = 1
      if (        eps(i,j,k) == 0 .and. eps(i+1,j,k) == 1 ) then
        tmp1(i+1,j,k) = C_05_R * tmp1(i+1,j,k)
        tmp2(i+1,j,k) = C_0_R
      end if
      do i = 2, imax-1
        if (      eps(i,j,k) == 0 .and. eps(i+1,j,k) == 1 ) then
          tmp1(i+1,j,k) = C_05_R * tmp1(i+1,j,k)
          tmp2(i+1,j,k) = C_0_R
        else if ( eps(i,j,k) == 0 .and. eps(i-1,j,k) == 1 ) then
          tmp1(i-1,j,k) = C_05_R * tmp1(i-1,j,k)
          tmp2(i-1,j,k) = C_0_R
        end if
      end do 
      i = imax
      if (        eps(i,j,k) == 0 .and. eps(i-1,j,k) == 1 ) then
        tmp1(i-1,j,k) = C_05_R * tmp1(i-1,j,k)
        tmp2(i-1,j,k) = C_0_R
      end if
    end do

    ! in y-direction - inner points - boundary points
    do i = 1, imax   
      j = 1        
      if (        eps(i,j,k) == 0 .and. eps(i,j+1,k) == 1 ) then
        tmp1(i,j+1,k) = C_05_R * tmp1(i,j+1,k)
        tmp2(i,j+1,k) = C_0_R
      end if 
      do j = 2, jmax-1
        if (      eps(i,j,k) == 0 .and. eps(i,j+1,k) == 1 ) then
          tmp1(i,j+1,k) = C_05_R * tmp1(i,j+1,k)
          tmp2(i,j+1,k) = C_0_R
        else if ( eps(i,j,k) == 0 .and. eps(i,j-1,k) == 1 ) then
          tmp1(i,j-1,k) = C_05_R * tmp1(i,j-1,k)
          tmp2(i,j-1,k) = C_0_R
        end if
      end do 
      j = jmax
      if (        eps(i,j,k) == 0 .and. eps(i,j-1,k) == 1 ) then
        tmp1(i,j-1,k) = C_05_R * tmp1(i,j-1,k)
        tmp2(i,j-1,k) = C_0_R
      end if
    end do  
  end do

  ! in z-direction - inner points - boundary points
  do i = 1, imax   
    do j = 1, jmax 
      k = 1        
      if (        eps(i,j,k) == 0 .and. eps(i,j,k+1) == 1 ) then
        tmp1(i,j,k+1) = C_05_R * tmp1(i,j,k+1)
        tmp2(i,j,k+1) = C_0_R
      end if 
      do k = 2, kmax-1
        if (      eps(i,j,k) == 0 .and. eps(i,j,k+1) == 1 ) then
          tmp1(i,j,k+1) = C_05_R * tmp1(i,j,k+1)
          tmp2(i,j,k+1) = C_0_R
        else if ( eps(i,j,k) == 0 .and. eps(i,j,k-1) == 1 ) then
          tmp1(i,j,k-1) = C_05_R * tmp1(i,j,k-1)
          tmp2(i,j,k-1) = C_0_R
        end if
      end do 
      k = kmax
      if (        eps(i,j,k) == 0 .and. eps(i,j,k-1) == 1 ) then
        tmp1(i,j,k-1) = C_05_R * tmp1(i,j,k-1)
        tmp2(i,j,k-1) = C_0_R
      end if
    end do
  end do

  tmp1(:,:,:) = (C_1_R - tmp1(:,:,:))

  ! horizontal average - compute gamma_f
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac, g(3)%jac, gamma_f, wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac, g(3)%jac, gamma_0, wrk1d, area)

  ! compute gamma_i, gamma_s
  do j = 1, jmax 
    if ( gamma_f(j) /= C_1_R ) then
      gamma_i(j) = C_1_R - gamma_0(j) - gamma_1(j)
      gamma_s(j) = C_1_R - gamma_f(j)
      res_1(j)   = C_1_R - gamma_s(j) - gamma_f(j)
      res_2(j)   = C_1_R - gamma_0(j) - gamma_i(j) - gamma_1(j)
    else
      gamma_i(j) = C_1_R;  gamma_s(j) = C_1_R
      res_1(j)   = C_0_R;  res_2(j)   = C_0_R
    end if
  end do

  ! constistency check of partitioning
  check_1 = sum(res_1); check_2 = sum(res_2) 
#ifdef USE_MPI
  dummy = check_1
  CALL MPI_ALLREDUCE(dummy, check_1, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  dummy = check_2
  CALL MPI_ALLREDUCE(dummy, check_2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
  if (check_1 > 1e-16 .or. check_2 > 1e-16) then
    call TLAB_WRITE_ASCII(efile, 'IBM_AVG_GAMMA. Sum of gammas is not unity.')
    call TLAB_STOP(DNS_ERROR_IBM_GAMMA)
  end if

  nullify(wrk1d, res_1, res_2)

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