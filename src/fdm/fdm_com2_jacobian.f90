!########################################################################
!# Compact FDMs for non-uniform grids from Lele, JCP, 1992, usign Jacobian
!#
!# Maximum stencil size according to the form
!#
!# f''_i + a_1(f''_{i-1}+f''_{i+1}) + a_2(f''_{i-2}+f''_{i+2})
!#   = b_1(f_{i+1}-2f_i+f_{i-1}) + b_2(f_{i+2}-2f_i+f_{i-2}) + b_3(f_{i+3}-2f_i+f_{i-3})
!#
!# biased at the boundary according to the form
!#
!# f''_1 + a_1 f''_2 + a_2 f''_3
!#   = b_1 f_1 + b_2 f_2 + b_3 f_3 + b_4 f_4
!#
!# a_1 f''_1 + f''_2 + a_2 f''_3 =
!#   = b_1 f_1 + b_2 f_2 + b_3 f_3 + b_4 f_4
!#
!# a_1 f''_2 + f''_3 + a_2 f''_4 =
!#   = b_1 f_1 + b_2 f_2 + b_3 f_3 + b_4 f_4+  b_5 f_5 + b_6 f_6
!#
!# The stencil in the rhs at the boundary can have 1 more point, which is stored in r11
!#
!# System normalized s.t. 1. off-diagonal in B is 1.
!#
!########################################################################
module FDM_Com2_Jacobian
    use TLAB_CONSTANTS, only: wp, wi, pi_wp
    implicit none
    private

    public :: FDM_C2N4_Jacobian         ! 2. order derivative, 4. order approximation
    public :: FDM_C2N6_Jacobian         ! 2. order derivative, 6. order approximation, hypodiffusive
    public :: FDM_C2N6_Hyper_Jacobian   ! 2. order derivative, 6. order approximation, hyperdiffusive

    logical periodic_loc

contains
    !########################################################################
    ! rhs is still pentadiagonal because of the bcs
    subroutine FDM_C2N4_Jacobian(nmax, dx, lhs, rhs, coef, periodic)
        integer(wi), intent(in) :: nmax
        real(wp), intent(in) :: dx(nmax, 2)
        real(wp), intent(out) :: lhs(nmax, 3)       ! LHS diagonals; a_2 = 0
        real(wp), intent(out) :: rhs(nmax, 5 + 3)   ! RHS diagonals; b_2, b_3 = 0
        real(wp), intent(out) :: coef(5)            ! a_1, a_2, b_1, b_2, b_3
        logical, intent(in), optional :: periodic

        ! -------------------------------------------------------------------
        real(wp) coef_bc1(6)

        if (present(periodic)) then
            periodic_loc = periodic
        else
            periodic_loc = .false.
        end if

        ! #######################################################################
        ! Interior points according to Eq. 2.2.7, 6th order approximation
        coef(1:2) = [0.1_wp, 0.0_wp]                        ! a_1, a_2
        coef(3:5) = [1.2_wp, 0.0_wp, 0.0_wp]                ! b_1, b_2, b_3

        if (periodic_loc) then
            call Create_System_2der(dx, lhs, rhs(:, 1:5), rhs(:, 6:8), coef)

        else    ! biased at the boundaries
            ! 3rd order, Eq. 4.3.1
            coef_bc1(1:2) = [11.0_wp, 0.0_wp]                       ! a_1, a_2
            coef_bc1(3:6) = [13.0_wp, -27.0_wp, 15.0_wp, -1.0_wp]   ! b_1, b_2, b_3, b_4

            call Create_System_2der(dx, lhs, rhs(:, 1:5), rhs(:, 6:8), coef, coef_bc1)

        end if

        return
    end subroutine FDM_C2N4_Jacobian

!########################################################################
    subroutine FDM_C2N6_Jacobian(nmax, dx, lhs, rhs, coef, periodic)
        integer(wi), intent(in) :: nmax
        real(wp), intent(in) :: dx(nmax, 2)
        real(wp), intent(out) :: lhs(nmax, 3)       ! LHS diagonals; a_2 = 0
        real(wp), intent(out) :: rhs(nmax, 5 + 3)   ! RHS diagonals; b_3 = 0
        real(wp), intent(out) :: coef(5)            ! a_1, a_2, b_1, b_2, b_3
        logical, intent(in), optional :: periodic

        ! -------------------------------------------------------------------
        real(wp) coef_bc1(6), coef_bc2(6)

        if (present(periodic)) then
            periodic_loc = periodic
        else
            periodic_loc = .false.
        end if

        ! #######################################################################
        ! Interior points according to Eq. 2.2.7, 6th order approximation
        coef(1:2) = [2.0_wp/11.0_wp, 0.0_wp]                        ! a_1, a_2
        coef(3:5) = [12.0_wp/11.0_wp, 3.0_wp/44.0_wp, 0.0_wp]       ! b_1, b_2, b_3

        if (periodic_loc) then
            call Create_System_2der(dx, lhs, rhs(:, 1:5), rhs(:, 6:8), coef)

        else    ! biased at the boundaries
            ! 3rd order, Eq. 4.3.1
            coef_bc1(1:2) = [11.0_wp, 0.0_wp]                       ! a_1, a_2
            coef_bc1(3:6) = [13.0_wp, -27.0_wp, 15.0_wp, -1.0_wp]   ! b_1, b_2, b_3, b_4

            ! 4th order, Eq. 2.2.6 (alpha=1/10).
            coef_bc2(1:2) = [0.1_wp, 0.1_wp]                        ! a_1, a_2
            coef_bc2(3:6) = [1.2_wp, -2.4_wp, 1.2_wp, 0.0_wp]       ! b_1, b_2, b_3, b_4

            call Create_System_2der(dx, lhs, rhs(:, 1:5), rhs(:, 6:8), coef, coef_bc1, coef_bc2)

        end if

        return
    end subroutine FDM_C2N6_Jacobian

!########################################################################
    subroutine FDM_C2N6_Hyper_Jacobian(nmax, dx, lhs, rhs, coef, periodic)
        integer(wi), intent(in) :: nmax
        real(wp), intent(in) :: dx(nmax, 2)
        real(wp), intent(out) :: lhs(nmax, 3)       ! LHS diagonals; a_2 = 0
        real(wp), intent(out) :: rhs(nmax, 7 + 3)   ! RHS diagonals; b_3 = 0
        real(wp), intent(out) :: coef(5)            ! a_1, a_2, b_1, b_2, b_3
        logical, intent(in), optional :: periodic

        ! -------------------------------------------------------------------
        real(wp) coef_bc1(6), coef_bc2(6), coef_bc3(8), kc

        if (present(periodic)) then
            periodic_loc = periodic
        else
            periodic_loc = .false.
        end if

        ! #######################################################################
        ! Interior points according to to JCP, Lamballais et al. 2011, JCP 230:3270-3275, Eqs. 1,3, 6th order
        ! One more diagonal in rhs to better match the transfer function (hyper- instead of hypodiffusive)
        kc = pi_wp**2.0_wp
        coef(1:2) = [(272.0_wp - 45.0_wp*kc)/(416.0_wp - 90.0_wp*kc), &
                     0.0_wp]                                                             ! a_1, a_2
        coef(3:5) = [(48.0_wp - 135.0_wp*kc)/(1664.0_wp - 360.0_wp*kc), &
                     (528.0_wp - 81.0_wp*kc)/(208.0_wp - 45.0_wp*kc)/4.0_wp, &
                     -(432.0_wp - 63.0_wp*kc)/(1664.0_wp - 360.0_wp*kc)/9.0_wp]     ! b_1, b_2, b_3

        if (periodic_loc) then
            call Create_System_2der(dx, lhs, rhs(:, 1:7), rhs(:, 8:10), coef)

        else    ! biased at the boundaries
            ! 3rd order, Eq. 4.3.1
            coef_bc1(1:2) = [11.0_wp, 0.0_wp]                       ! a_1, a_2
            coef_bc1(3:6) = [13.0_wp, -27.0_wp, 15.0_wp, -1.0_wp]   ! b_1, b_2, b_3, b_4

            ! 4th order, Eq. 2.2.6 (alpha=1/10).
            coef_bc2(1:2) = [0.1_wp, 0.1_wp]                        ! a_1, a_2
            coef_bc2(3:6) = [1.2_wp, -2.4_wp, 1.2_wp, 0.0_wp]       ! b_1, b_2, b_3, b_4

            ! 6th order, Eq. 2.2.7, 6th order approximation
            coef_bc3(1:2) = [2.0_wp/11.0_wp, 2.0_wp/11.0_wp]                ! a_1, a_2
            coef_bc3(3:8) = [3.0_wp/44.0_wp, 12.0_wp/11.0_wp, -51.0_wp/22.0_wp, 12.0_wp/11.0_wp, 3.0_wp/44.0_wp, 0.0_wp]       ! b_1, b_2, b_3, b_4, b_5, b_6

            call Create_System_2der(dx, lhs, rhs(:, 1:7), rhs(:, 8:10), coef, coef_bc1, coef_bc2, coef_bc3)

        end if

        return
    end subroutine FDM_C2N6_Hyper_Jacobian

!########################################################################
    subroutine Create_System_2der(dx, lhs, rhs, rhs_d1, coef_int, coef_bc1, coef_bc2, coef_bc3)
        real(wp), intent(in) :: dx(:, :)        ! 1. and 2. order Jacobians
        real(wp), intent(out) :: lhs(:, :)      ! LHS diagonals
        real(wp), intent(out) :: rhs(:, :)      ! RHS diagonals
        real(wp), intent(out) :: rhs_d1(:, :)   ! RHS diagonals that multiply 1. derivative
        real(wp), intent(in) :: coef_int(5)             ! a_1, a_2b_1, b_2, b_3
        real(wp), intent(in), optional :: coef_bc1(6)   ! a_1, a_2, b_1, b_2, b_3, b_4
        real(wp), intent(in), optional :: coef_bc2(6)   ! a_1, a_2, b_1, b_2, b_3, b_4
        real(wp), intent(in), optional :: coef_bc3(8)   ! a_1, a_2, b_1, b_2, b_3, b_4, b_5, b_6

        ! -------------------------------------------------------------------
        integer(wi) n, nmax, idl, idr, i, imax

        ! #######################################################################
        idl = size(lhs, 2)/2 + 1        ! center diagonal in lhs
        idr = size(rhs, 2)/2 + 1        ! center diagonal in rhs
        nmax = size(lhs, 1)             ! # grid points

        ! lhs diagonals
        lhs(:, idl) = 1.0_wp                ! center diagonal
        do i = 1, idl - 1                   ! off-diagonals
            lhs(:, idl - i) = coef_int(i)
            lhs(:, idl + i) = coef_int(i)
        end do

        ! lhs diagonals
        rhs(:, idr) = 0.0_wp                ! initialize center diagonal
        do i = 1, idr - 1
            rhs(:, idr) = rhs(:, idr) - 2.0_wp*coef_int(i + 2)  ! center diagonal
            rhs(:, idr - i) = coef_int(i + 2)                   ! off-diagonals
            rhs(:, idr + i) = coef_int(i + 2)
        end do

        ! boundaries
        if (present(coef_bc1)) then
            n = 1
            lhs(n, idl) = 1.0_wp                                    ! lhs center diagonal
            if (idl > 1) then
                imax = min(idl - 1, 2)                              ! max of 3 point stencil, the first one set to 1
                lhs(n, idl + 1:idl + imax) = coef_bc1(1:imax)       ! lhs off-diagonals
            end if
            rhs(n, :) = 0.0_wp
            imax = min(idr, 4)                                      ! max of 4 point stencil
            rhs(n, idr:idr + imax - 1) = coef_bc1(3:3 + imax - 1)   ! rhs center and off-diagonals
            rhs(n, 1) = coef_bc1(3 + imax)                          ! extended rhs stencil

            n = nmax        ! symmetry property to define values at end
            lhs(n, :) = lhs(1, size(lhs, 2):1:-1)
            rhs(n, :) = rhs(1, size(rhs, 2):1:-1)
        end if

        if (present(coef_bc2)) then
            n = 2
            if (size(lhs, 2) == 3) then
                lhs(n, :) = [coef_bc2(1), 1.0_wp, coef_bc2(2)]
            else
                lhs(n, :) = [1.0_wp]
            end if
            rhs(n, :) = 0.0_wp
            imax = min(idr + 1, 4)                                      ! max of 4 point stencil
            rhs(n, idr - 1:idr + imax - 2) = coef_bc2(3:3 + imax - 1)   ! rhs center and off-diagonals

            n = nmax - 1    ! symmetry property to define values at end
            lhs(n, :) = lhs(2, size(lhs, 2):1:-1)
            rhs(n, :) = rhs(2, size(rhs, 2):1:-1)
        end if

        if (present(coef_bc3)) then
            n = 3
            if (size(lhs, 2) == 3) then
                lhs(n, :) = [coef_bc3(1), 1.0_wp, coef_bc3(2)]
            else
                lhs(n, :) = [1.0_wp]
            end if
            rhs(n, :) = 0.0_wp
            imax = min(idr + 2, 6)                                      ! max of 6 point stencil
            rhs(n, idr - 2:idr + imax - 3) = coef_bc3(3:3 + imax - 1)   ! rhs center and off-diagonals

            n = nmax - 2    ! symmetry property to define values at end
            lhs(n, :) = lhs(3, size(lhs, 2):1:-1)
            rhs(n, :) = rhs(3, size(rhs, 2):1:-1)
        end if

        ! multiply by the Jacobians
        rhs_d1(:, idl) = -lhs(:, idl)*dx(:, 2)          ! center diagonal
        do i = 1, idl - 1                               ! off-diagonals
            rhs_d1(:, idl - i) = -lhs(:, idl - i)*cshift(dx(:, 2), -i)
            rhs_d1(:, idl + i) = -lhs(:, idl + i)*cshift(dx(:, 2), +i)
        end do

        lhs(:, idl) = lhs(:, idl)*dx(:, 1)*dx(:, 1)     ! center diagonal
        do i = 1, idl - 1                               ! off-diagonals
            lhs(:, idl - i) = lhs(:, idl - i)*cshift(dx(:, 1), -i)*cshift(dx(:, 1), -i)
            lhs(:, idl + i) = lhs(:, idl + i)*cshift(dx(:, 1), +i)*cshift(dx(:, 1), +i)
        end do

        ! normalize such the coefficent in 1. off-diagonal of rhs is 1
        lhs(:, :) = lhs(:, :)/coef_int(3)
        rhs(:, :) = rhs(:, :)/coef_int(3)
        rhs_d1(:, :) = rhs_d1(:, :)/coef_int(3)

        return
    end subroutine Create_System_2der

end module FDM_Com2_Jacobian
