!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/16 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Coefficients of pentadiagonal 6th-order scheme for 1st derivative.
!# WARNING: duplicated code, cf. coeff in FDM_C1N6_Jacobian_Penta!
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################

subroutine FDM_C1N6M_COEFF()
    use TLab_Constants, only: wp, wi

    use FDM_PROCS, only: C1N6M_ALPHA, C1N6M_BETA
    use FDM_PROCS, only: C1N6M_ALPHA2, C1N6M_BETA2
    use FDM_PROCS, only: C1N6M_A, C1N6M_B, C1N6M_C
    use FDM_PROCS, only: C1N6M_AD2, C1N6M_BD4, C1N6M_CD6
    use FDM_PROCS, only: C1N6M_BD2, C1N6M_CD3

    implicit none

! #######################################################################
! -------------------------------------------------------------------
! Pentadiagonal 6th-order scheme (JCP Lele 1992, Eq. 2.1.10) with
! similar truncation error as 6th-order tridiagonal scheme
! (JCP Lele 1992, Eq. 2.1.7 with alpha=1/3). Value of beta is derived
! with table 1 (2.1.7 = 2.1.10). Largest stable value with alpha=0.56.
! (Simulations are unstable for larger alpha values!)
! (values of alpha and beta can be changed, according to Lele...)
! -------------------------------------------------------------------
    C1N6M_ALPHA = 56.0_wp/100.0_wp
    C1N6M_BETA = (2.0_wp/5.0_wp)*(-1.0_wp/3.0_wp + C1N6M_ALPHA)
    C1N6M_A = (1.0_wp/6.0_wp)*(9.0_wp + C1N6M_ALPHA - 20.0_wp*C1N6M_BETA)
    C1N6M_B = (1.0_wp/15.0_wp)*(-9.0_wp + 32.0_wp*C1N6M_ALPHA + 62.0_wp*C1N6M_BETA)
    C1N6M_C = (1.0_wp/10.0_wp)*(1.0_wp - 3.0_wp*C1N6M_ALPHA + 12.0_wp*C1N6M_BETA)

! #######################################################################
    C1N6M_AD2 = C1N6M_A/2.0_wp
    C1N6M_BD4 = C1N6M_B/4.0_wp
    C1N6M_CD6 = C1N6M_C/6.0_wp
    !
    C1N6M_BD2 = C1N6M_B/2.0_wp
    C1N6M_CD3 = C1N6M_C/3.0_wp
    !
    C1N6M_ALPHA2 = 2.0_wp*C1N6M_ALPHA
    C1N6M_BETA2 = 2.0_wp*C1N6M_BETA

    return
end subroutine FDM_C1N6M_COEFF
