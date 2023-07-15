!########################################################################
!# Compact FDMs for non-uniform grids from Lele, JCP, 1992, usign Jacobian
!#
!# Maximum stencil size according to the form
!#
!# f'_i + a_1(f'_{i-1}+f'_{i+1}) + a_2(f'_{i-2}+f'_{i+2})
!#   = b_1(f_{i+1}-f_{i-1}) + b_2(f_{i+2}-f_{i-2}) + b_3(f_{i+3}-f_{i-3})
!#
!# biased at the boundary according to the form
!#
!# f'_1 + a_1 f'_2 + a_2 f'_3
!#   = b_1 f_1 + b_2 f_2 + b_3 f_3 + b_4 f_4
!#
!# a_1 f'_1 + f'_2 + a_2 f'_3 =
!#   = b_1 f_1 + b_2 f_2 + b_3 f_3 + b_4 f_4
!#
!# The stencil in the rhs at the boundary can have 1 more point, which is stored in r11
!#
!# System normalized s.t. 1. off-diagonal in B is 1.
!#
!########################################################################
module FDM_Com1_Jacobian
    use TLAB_CONSTANTS, only: wp, wi
    implicit none
    private

    public :: FDM_C1N4_Jacobian     ! 1. order derivative, 4. order approximation
    public :: FDM_C1N6_Jacobian     ! 1. order derivative, 6. order approximation

    logical periodic_loc

contains
    !########################################################################
    subroutine FDM_C1N4_Jacobian(nx, dx, lhs, rhs, coef, periodic)
        integer(wi), intent(in) :: nx
        real(wp), intent(in) :: dx(nx)
        real(wp), intent(out) :: lhs(nx, 3)     ! LHS diagonals; a_2 = 0
        real(wp), intent(out) :: rhs(nx, 3)     ! RHS diagonals; b_2, b_3 = 0
        real(wp), intent(out) :: coef(5)        ! a_1, a_2, b_1, b_2, b_3
        logical, intent(in), optional :: periodic

        ! -------------------------------------------------------------------
        real(wp) coef_bc1(6)

        if (present(periodic)) then
            periodic_loc = periodic
        else
            periodic_loc = .false.
        end if

        ! #######################################################################
        ! Interior points according to Eq. 2.1.7 (\alpha=1/3), 6th order approximation
        coef(1:2) = [0.25_wp, 0.0_wp]                                   ! a_1, a_2
        coef(3:5) = [0.75_wp, 0.0_wp, 0.0_wp]                           ! b_1, b_2, b_3

        if (periodic_loc) then
            call Create_System_1der(dx, lhs, rhs, coef)

        else    ! biased at the boundaries
            ! 3rd order, Eq. 4.1.3 with \alpha=2
            coef_bc1(1:2) = [2.0_wp, 0.0_wp]                            ! a_1, a_2
            coef_bc1(3:6) = [-15.0_wp/6.0_wp, 2.0_wp, 0.5_wp, 0.0_wp]   ! b_1, b_2, b_3, b_4

            call Create_System_1der(dx, lhs, rhs, coef, coef_bc1)

        end if

        return
    end subroutine FDM_C1N4_Jacobian

!########################################################################
    subroutine FDM_C1N6_Jacobian(nx, dx, lhs, rhs, coef, periodic)
        integer(wi), intent(in) :: nx
        real(wp), intent(in) :: dx(nx)
        real(wp), intent(out) :: lhs(nx, 3)     ! LHS diagonals; a_2 = 0
        real(wp), intent(out) :: rhs(nx, 5)     ! RHS diagonals; b_3 = 0
        real(wp), intent(out) :: coef(5)        ! a_1, a_2, b_1, b_2, b_3
        logical, intent(in), optional :: periodic

        ! -------------------------------------------------------------------
        real(wp) coef_bc1(6), coef_bc2(6)

        if (present(periodic)) then
            periodic_loc = periodic
        else
            periodic_loc = .false.
        end if

        ! #######################################################################
        ! Interior points according to Eq. 2.1.7 (\alpha=1/3), 6th order approximation
        coef(1:2) = [1.0_wp/3.0_wp, 0.0_wp]                                     ! a_1, a_2
        coef(3:5) = [7.0_wp/9.0_wp, 1.0_wp/36.0_wp, 0.0_wp]                     ! b_1, b_2, b_3

        if (periodic_loc) then
            call Create_System_1der(dx, lhs, rhs, coef)

        else    ! biased at the boundaries
            ! 3rd order, Eq. 4.1.3 with \alpha=2
            coef_bc1(1:2) = [2.0_wp, 0.0_wp]                                    ! a_1, a_2
            coef_bc1(3:6) = [-15.0_wp/6.0_wp, 2.0_wp, 0.5_wp, 0.0_wp]           ! b_1, b_2, b_3, b_4

            ! 5th order in Carpenter et al, JCP, 108:272-295, 1993, Eq. 93, who study the effect of boundary points on stability
            coef_bc2(1:2) = [1.0_wp/6.0_wp, 0.5_wp]                             ! a_1, a_2
            coef_bc2(3:6) = [-5.0_wp/9.0_wp, -0.5_wp, 1.0_wp, 1.0_wp/18.0_wp]   ! b_1, b_2, b_3, b_4
            ! 4th order, Eq. 2.1.6 with \alpha=1/4.
            ! coef_bc2(1:2) = [0.25_wp, 0.25_wp]                                ! a_1, a_2
            ! coef_bc2(3:6) = [-0.75_wp, 0.0_wp, 0.75_wp, 0.0_wp]               ! b_1, b_2, b_3, b_4

            call Create_System_1der(dx, lhs, rhs, coef, coef_bc1, coef_bc2)

        end if

        return
    end subroutine FDM_C1N6_Jacobian

!########################################################################
    subroutine Create_System_1der(dx, lhs, rhs, coef_int, coef_bc1, coef_bc2)
        real(wp), intent(in) :: dx(:)
        real(wp), intent(out) :: lhs(:, :)   ! LHS diagonals
        real(wp), intent(out) :: rhs(:, :)   ! RHS diagonals
        real(wp), intent(in) :: coef_int(5)             ! a_1, a_2b_1, b_2, b_3
        real(wp), intent(in), optional :: coef_bc1(6)   ! a_1, a_2, b_1, b_2, b_3, b_4
        real(wp), intent(in), optional :: coef_bc2(6)   ! a_1, a_2, b_1, b_2, b_3, b_4

        ! -------------------------------------------------------------------
        integer(wi) n, nx, idl, idr, i, imax

        ! #######################################################################
        idl = size(lhs, 2)/2 + 1        ! center diagonal in lhs
        idr = size(rhs, 2)/2 + 1        ! center diagonal in rhs
        nx = size(lhs, 1)             ! # grid points

        ! lhs diagonals
        lhs(:, idl) = 1.0_wp                ! center diagonal
        do i = 1, idl - 1                   ! off-diagonals
            lhs(:, idl - i) = coef_int(i)
            lhs(:, idl + i) = coef_int(i)
        end do

        ! lhs diagonals
        rhs(:, idr) = 0.0_wp                ! center diagonal
        do i = 1, idr - 1                   ! off-diagonals
            rhs(:, idr - i) = -coef_int(i + 2)
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

            n = nx                                                  ! symmetry property to define values at end
            lhs(n, :) = lhs(1, size(lhs, 2):1:-1)
            rhs(n, :) = -rhs(1, size(rhs, 2):1:-1)
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

            n = nx - 1                                                  ! symmetry property to define values at end
            lhs(n, :) = lhs(2, size(lhs, 2):1:-1)
            rhs(n, :) = -rhs(2, size(rhs, 2):1:-1)
        end if

        ! multiply by the Jacobian
        lhs(:, idl) = lhs(:, idl)*dx(:)     ! center diagonal
        do i = 1, idl - 1                   ! off-diagonals
            lhs(:, idl - i) = lhs(:, idl - i)*cshift(dx(:), -i)
            lhs(:, idl + i) = lhs(:, idl + i)*cshift(dx(:), +i)
        end do

        ! normalize such the coefficent in 1. off-diagonal of rhs is 1
        lhs(:, :) = lhs(:, :)/coef_int(3)
        rhs(:, :) = rhs(:, :)/coef_int(3)

        return
    end subroutine Create_System_1der

end module FDM_Com1_Jacobian
