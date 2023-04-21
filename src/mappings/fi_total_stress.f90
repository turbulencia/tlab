#include "dns_const.h"
!########################################################################
!# HISTORY / AUTHORS
!#
!# 2023/03/28 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   Calculate the total stress tensor
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

module FI_TOTAL_STRESS
    
    use TLAB_CONSTANTS, only: wp, wi
    use TLAB_VARS,      only: g, visc, imode_ibm
    use IBM_VARS,       only: ibm_partial
    use OPR_PARTIAL

    implicit none
    
    private

    integer(wi), parameter :: bcs(2, 2) = 0

    public :: FI_TOTAL_STRESS_TENSOR

contains

!########################################################################

    subroutine FI_TOTAL_STRESS_TENSOR(nx, ny, nz, u, v, w, p, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)
        
        integer(wi),                   intent(in   ) :: nx, ny, nz
        real(wp), dimension(nx,ny,nz), intent(in   ) :: u, v, w, p
        real(wp), dimension(nx,ny,nz), intent(  out) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 ! Sxx,Syy,Szz,Sxy,Sxz,Syz
        
        ! ================================================================== !
    
        ! OPR_PARTIAL_X/Y/Z uses modified fields for derivatives
        if ( imode_ibm == 1 ) ibm_partial = .true.

        ! Off diagonal terms

        ! Uy, Vx
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), v, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp4)
        tmp4 = visc*(tmp4 + tmp1)

        ! Uz, Wx
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp1)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, tmp5)
        tmp5 = visc*(tmp5 + tmp1)

        ! Vz, Wy
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp1)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp6)
        tmp6 = visc*(tmp6 + tmp1)

        ! diagonal terms
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), u, tmp1)
        tmp1 =  -p !+ 2.0_wp*visc*tmp1
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), v, tmp2)
        tmp2 =  -p !+ 2.0_wp*visc*tmp2
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), w, tmp3)
        tmp3 =  -p !+ 2.0_wp*visc*tmp3

        if ( imode_ibm == 1 ) ibm_partial = .false.

        return

    end subroutine FI_TOTAL_STRESS_TENSOR
    
end module FI_TOTAL_STRESS