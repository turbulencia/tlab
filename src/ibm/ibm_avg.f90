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
!#   ==> gamma_s & gamma_f are computed afterwards in python!
!#       (not implemented in fortran, challenging with MPI ...) 
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

subroutine IBM_AVG_GAMMA(gamma_0, gamma_1, eps, tmp1)
    use TLab_Constants, only: wp
    use TLab_Memory, only: imax, jmax, kmax, isize_field
    use Averages, only: AVG_IK_V

    implicit none

    real(wp), dimension(jmax), intent(out) :: gamma_0, gamma_1
    real(wp), dimension(imax, jmax, kmax), intent(in) :: eps
    real(wp), dimension(isize_field), intent(inout) :: tmp1

    ! ================================================================== !
    ! horizontal average - compute gamma_1
    call AVG_IK_V(imax, jmax, kmax, eps, gamma_1, tmp1(1:jmax))

    gamma_0(:) = 1.0_wp - gamma_1(:) 

    return
end subroutine IBM_AVG_GAMMA
!########################################################################

subroutine IBM_AVG_SCAL_BCS(is, scalv_bcs)

    use IBM_VARS, only: ibm_objup, max_height_objlo, max_height_objup
    use IBM_VARS, only: ibmscaljmin, ibmscaljmax
    use TLab_Memory, only: jmax
    use TLab_Constants, only: wi, wp

    implicit none

    integer(wi), intent(in) :: is
    real(wp), dimension(jmax), intent(out) :: scalv_bcs

    integer(wi) :: j

    ! ================================================================== !
    ! write out scalar boundary values applied in solids (vertical)
    ! assuming homogenous temperature in solids, otherwise change here

    scalv_bcs = 0.0_wp

    do j = 1, int(max_height_objlo, wi)
        scalv_bcs(j) = ibmscaljmin(is)
    end do
    if (ibm_objup) then
        do j = jmax - int(max_height_objup, wi), jmax
            scalv_bcs(j) = ibmscaljmax(is)
        end do
    end if

    return
end subroutine IBM_AVG_SCAL_BCS
