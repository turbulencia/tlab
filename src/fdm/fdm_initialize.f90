#include "dns_const.h"

!# Compute dx/di and create LU factorization for first- and second-order derivatives

subroutine FDM_INITIALIZE(x, g, wrk1d)
    use TLAB_CONSTANTS, only: wp, wi, pi_wp
#ifdef TRACE_ON
    use TLAB_CONSTANTS, only: tfile
#endif
    use TLAB_TYPES, only: grid_dt
    use TLAB_VARS, only: inb_scal, stagger_on
    use TLAB_VARS, only: reynolds, schmidt
    use TLAB_VARS, only: C1N6M_ALPHA2, C1N6M_BETA2
    use TLAB_VARS, only: C1N6M_A, C1N6M_BD2, C1N6M_CD3
    use TLAB_PROCS
    use FDM_COM_DIRECT

    implicit none

    type(grid_dt), intent(inout) :: g
    real(wp), intent(inout) :: x(g%size, g%inb_grid)
    real(wp), intent(inout) :: wrk1d(g%size, 5)

    target x

! -------------------------------------------------------------------
    integer(wi) i, ip, is, ig, ibc_min, ibc_max, nx
    real(wp) dummy
    real(wp) a1, a2, b1, b2, b3, kc

    integer, parameter :: i0 = 0, i1 = 1

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'Entering SUBROUTINE FDM_INITIALIZE')
#endif

! ###################################################################
    nx = g%size

    ig = 1 ! Accumulating counter to define pointers inside array x

! ###################################################################
    g%nodes => x(:, ig)

    ig = ig + 1

! ###################################################################
! Coefficients of pentadiagonal 6th-order scheme for 1st derivative
    if (g%mode_fdm == FDM_COM6_JACPENTA) call FDM_C1N6M_COEFF()

! ###################################################################
! Jacobians
! ###################################################################
    g%jac => x(:, ig:)

    if (nx == 1) then
        g%jac(1, 1) = 1.0_wp
        g%jac(1, 2) = 1.0_wp
        g%jac(1, 3) = 0.0_wp
        g%jac(1, 4) = 0.0_wp
        return
    end if

! ###################################################################
    if (g%uniform) then
! -------------------------------------------------------------------
! first derivative
! -------------------------------------------------------------------
        do i = 2, nx - 1
            g%jac(i, 1) = (x(i + 1, 1) - x(i - 1, 1))*0.5_wp
        end do

! Boundary points
        if (g%periodic) then
            g%jac(nx, 1) = (x(1, 1) + g%scale - x(nx - 1, 1))*0.5_wp
            g%jac(1, 1) = (x(2, 1) - x(1, 1) + x(1, 1) + g%scale - x(nx, 1))*0.5_wp

        else
            g%jac(1, 1) = g%jac(2, 1)
            g%jac(nx, 1) = g%jac(nx - 1, 1)

        end if

! -------------------------------------------------------------------
! second derivative is zero
! -------------------------------------------------------------------
        g%jac(:, 2) = 0.0_wp

! ###################################################################
    else ! derivative wrt computational grid, uniform
! -------------------------------------------------------------------
! first derivative
! -------------------------------------------------------------------
        g%jac(:, 1) = 1.0_wp

        select case (g%mode_fdm)

        case (FDM_COM4_JACOBIAN)
            call FDM_C1N4_LHS(nx, i0, i0, g%jac, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call FDM_C1N4_RHS(nx, i1, i0, i0, x, g%jac(1, 1))

        case (FDM_COM6_JACOBIAN, FDM_COM6_DIRECT)
            call FDM_C1N6_LHS(nx, i0, i0, g%jac, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call FDM_C1N6_RHS(nx, i1, i0, i0, x, g%jac(1, 1))

        case (FDM_COM6_JACPENTA)
            call FDM_C1N6M_LHS(nx, i0, i0, g%jac, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
            call FDM_C1N6M_RHS(nx, i1, i0, i0, x, g%jac(1, 1))

        case (FDM_COM8_JACOBIAN)
            call FDM_C1N8_LHS(nx, i0, i0, g%jac, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call FDM_C1N8_RHS(nx, i1, i0, i0, x, g%jac(1, 1))

        end select

        if (.not. (g%mode_fdm == FDM_COM6_JACPENTA)) then
            call TRIDFS(nx, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call TRIDSS(nx, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), g%jac(1, 1))
        else
            call PENTADFS2(nx, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
            call PENTADSS2(nx, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), g%jac(1, 1))
        end if

! -------------------------------------------------------------------
! second derivative
! -------------------------------------------------------------------
        wrk1d(:, 4) = 1.0_wp; wrk1d(:, 5) = 0.0_wp

        select case (g%mode_fdm)

        case (FDM_COM4_JACOBIAN)
            call FDM_C2N4_LHS(nx, i0, i0, wrk1d(1, 4), wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call FDM_C2N4_RHS(nx, i1, i0, i0, x, g%jac(1, 2))

        case (FDM_COM6_JACOBIAN, FDM_COM6_DIRECT, FDM_COM6_JACPENTA)
            ! call FDM_C2N6_LHS(nx, i0, i0, wrk1d(1, 4), wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            ! call FDM_C2N6_RHS(nx, i1, i0, i0, x, g%jac(1, 2))
            call FDM_C2N6H_LHS(nx, i0, i0, wrk1d(1, 4), wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call FDM_C2N6H_RHS(nx, i1, i0, i0, x, g%jac(1, 2))

        case (FDM_COM8_JACOBIAN) ! Not yet developed; default to 6. order
            call FDM_C2N6_LHS(nx, i0, i0, wrk1d(1, 4), wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call FDM_C2N6_RHS(nx, i1, i0, i0, x, g%jac(1, 2))

        end select

        call TRIDFS(nx, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
        call TRIDSS(nx, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), g%jac(1, 2))

    end if

! ###################################################################
! Saving operations for the time-stability constraint
    g%jac(:, 3) = 1.0_wp/g%jac(:, 1)
    g%jac(:, 4) = g%jac(:, 3)*g%jac(:, 3)

    ig = ig + 4

! ###################################################################
! LU factorization first-order derivative, done in routine TRID*FS
! ###################################################################
    g%lu1 => x(:, ig:)

! -------------------------------------------------------------------
! Periodic case
! -------------------------------------------------------------------
    if (g%periodic) then
        select case (g%mode_fdm)

        case (FDM_COM4_JACOBIAN)
            call FDM_C1N4P_LHS(nx, g%jac, g%lu1(1, 1), g%lu1(1, 2), g%lu1(1, 3))

        case (FDM_COM6_JACOBIAN, FDM_COM6_DIRECT) ! Direct = Jacobian because uniform grid
            call FDM_C1N6P_LHS(nx, g%jac, g%lu1(1, 1), g%lu1(1, 2), g%lu1(1, 3))

        case (FDM_COM6_JACPENTA)
            call FDM_C1N6MP_LHS(nx, g%jac, g%lu1(1, 1), g%lu1(1, 2), g%lu1(1, 3), g%lu1(1, 4), g%lu1(1, 5))

        case (FDM_COM8_JACOBIAN)
            call FDM_C1N8P_LHS(nx, g%jac, g%lu1(1, 1), g%lu1(1, 2), g%lu1(1, 3))

        end select

        if (.not. (g%mode_fdm == FDM_COM6_JACPENTA)) then
            call TRIDPFS(nx, g%lu1(1, 1), g%lu1(1, 2), g%lu1(1, 3), g%lu1(1, 4), g%lu1(1, 5))
        else
            call PENTADPFS(nx, g%lu1(1, 1), g%lu1(1, 2), g%lu1(1, 3), g%lu1(1, 4), g%lu1(1, 5), g%lu1(1, 6), g%lu1(1, 7))
        end if
        ig = ig + 7

! -------------------------------------------------------------------
! Nonperiodic case (4 different BCs)
! -------------------------------------------------------------------
    else
        do i = 0, 3
            ibc_min = mod(i, 2)
            ibc_max = i/2
            ip = i*5
            select case (g%mode_fdm)

            case (FDM_COM4_JACOBIAN)
                call FDM_C1N4_LHS(nx, ibc_min, ibc_max, g%jac, g%lu1(1, ip + 1), g%lu1(1, ip + 2), g%lu1(1, ip + 3))

            case (FDM_COM6_JACOBIAN)
                call FDM_C1N6_LHS(nx, ibc_min, ibc_max, g%jac, g%lu1(1, ip + 1), g%lu1(1, ip + 2), g%lu1(1, ip + 3))

            case (FDM_COM6_JACPENTA)
                call FDM_C1N6M_LHS(nx, ibc_min, ibc_max, g%jac, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3),g%lu1(1,ip+4),g%lu1(1,ip+5))

            case (FDM_COM8_JACOBIAN)
                call FDM_C1N8_LHS(nx, ibc_min, ibc_max, g%jac, g%lu1(1, ip + 1), g%lu1(1, ip + 2), g%lu1(1, ip + 3))

            case (FDM_COM6_DIRECT) ! Not yet implemented; using Jacobian version
                call FDM_C1N6_LHS(nx, ibc_min, ibc_max, g%jac, g%lu1(1, ip + 1), g%lu1(1, ip + 2), g%lu1(1, ip + 3))

            end select

            if (.not. (g%mode_fdm == FDM_COM6_JACPENTA)) then
                call TRIDFS(nx, g%lu1(1, ip + 1), g%lu1(1, ip + 2), g%lu1(1, ip + 3))
            else
                call PENTADFS2(nx, g%lu1(1, ip + 1), g%lu1(1, ip + 2), g%lu1(1, ip + 3), g%lu1(1, ip + 4), g%lu1(1, ip + 5))
            end if
            ig = ig + 5
        end do

    end if

! ###################################################################
! LU factorization second-order derivative, done in routine TRID*FS
! ###################################################################
    g%lu2 => x(:, ig:)

! -------------------------------------------------------------------
! Periodic case; pentadiagonal
! -------------------------------------------------------------------
    if (g%periodic) then
        select case (g%mode_fdm)

        case (FDM_COM4_JACOBIAN)
            call FDM_C2N4P_LHS(nx, g%jac, g%lu2(1, 1), g%lu2(1, 2), g%lu2(1, 3))

        case (FDM_COM6_JACOBIAN, FDM_COM6_DIRECT, FDM_COM6_JACPENTA)  ! Direct = Jacobian because uniform grid
            ! call FDM_C2N6P_LHS(nx, g%jac, g%lu2(1, 1), g%lu2(1, 2), g%lu2(1, 3))
            call FDM_C2N6HP_LHS(nx, g%jac, g%lu2(1, 1), g%lu2(1, 2), g%lu2(1, 3))

        case (FDM_COM8_JACOBIAN)                   ! Not yet developed
            call FDM_C2N6P_LHS(nx, g%jac, g%lu2(1, 1), g%lu2(1, 2), g%lu2(1, 3))

        end select

        call TRIDPFS(nx, g%lu2(1, 1), g%lu2(1, 2), g%lu2(1, 3), g%lu2(1, 4), g%lu2(1, 5))
        ig = ig + 5

! -------------------------------------------------------------------
! Nonperiodic case; tridiagonal for 4 different BCs
! -------------------------------------------------------------------
    else
        do i = 0, 3
            ibc_min = mod(i, 2)
            ibc_max = i/2
            ip = i*3
            select case (g%mode_fdm)

            case (FDM_COM4_JACOBIAN)
                call FDM_C2N4_LHS(nx, ibc_min, ibc_max, g%jac, g%lu2(1, ip + 1), g%lu2(1, ip + 2), g%lu2(1, ip + 3))

            case (FDM_COM6_JACOBIAN, FDM_COM6_JACPENTA)
                ! call FDM_C2N6_LHS(nx, ibc_min, ibc_max, g%jac, g%lu2(1, ip + 1), g%lu2(1, ip + 2), g%lu2(1, ip + 3))
                call FDM_C2N6H_LHS(nx, ibc_min, ibc_max, g%jac, g%lu2(1, ip + 1), g%lu2(1, ip + 2), g%lu2(1, ip + 3))

            case (FDM_COM8_JACOBIAN) ! Not yet implemented
                call FDM_C2N6_LHS(nx, ibc_min, ibc_max, g%jac, g%lu2(1, ip + 1), g%lu2(1, ip + 2), g%lu2(1, ip + 3)) ! 8th not yet developed

            case (FDM_COM6_DIRECT)
                if (i == 0) call FDM_C2N6ND_INITIALIZE(nx, x, g%lu2(:, ip + 1), g%lu2(:, ip + 4))

            case (FDM_COM4_DIRECT)
                if (i == 0) call FDM_C2N4ND_INITIALIZE(nx, x, g%lu2(:, ip + 1), g%lu2(:, ip + 4))

            end select

! The direct mode is only implemented for bcs=(0,0); we use the remaining array to save other data
            if (any([FDM_COM4_DIRECT, FDM_COM6_DIRECT] == g%mode_fdm)) then
                if (i == 0) then
                    g%lu2(:, ip + 8:ip + 10) = g%lu2(:, ip + 1:ip + 3) ! saving the array A w/o LU decomposition
                    call TRIDFS(nx, g%lu2(1, ip + 1), g%lu2(1, ip + 2), g%lu2(1, ip + 3))
                    ig = ig + 10
                end if
            else
                call TRIDFS(nx, g%lu2(1, ip + 1), g%lu2(1, ip + 2), g%lu2(1, ip + 3))
                ig = ig + 3
            end if

        end do

    end if

! ###################################################################
! LU factorization second-order derivative times the diffusivities
! ###################################################################
    g%lu2d => x(:, ig:)

    ip = 0
    do is = 0, inb_scal ! case 0 for the reynolds number
        if (is == 0) then; dummy = reynolds
        else; dummy = reynolds*schmidt(is); end if

! -------------------------------------------------------------------
! Periodic case; pentadiagonal
! -------------------------------------------------------------------
        if (g%periodic) then ! Check routines TRIDPFS and TRIDPSS
            g%lu2d(:, ip + 1) = g%lu2(:, 1)             ! matrix L; 1. subdiagonal
            g%lu2d(:, ip + 2) = g%lu2(:, 2)/dummy       ! matrix L; 1/diagonal
            g%lu2d(:, ip + 3) = g%lu2(:, 3)             ! matrix U is the same
            g%lu2d(:, ip + 4) = g%lu2(:, 4)*dummy       ! matrix L; Additional row/column
            g%lu2d(:, ip + 5) = g%lu2(:, 5)             ! matrix U is the same

            ig = ig + 5
            ip = ip + 5

! -------------------------------------------------------------------
! Nonperiodic case; tridiagonal, 1 single BCs
! -------------------------------------------------------------------
        else                   ! Check routines TRIDFS and TRIDSS
            g%lu2d(:, ip + 1) = g%lu2(:, 1)             ! matrix L is the same
            g%lu2d(:, ip + 2) = g%lu2(:, 2)/dummy       ! matrix U; 1/diagonal
            g%lu2d(:, ip + 3) = g%lu2(:, 3)*dummy       ! matrix U; 1. superdiagonal

            ig = ig + 3
            ip = ip + 3

        end if

    end do

! ###################################################################
! LU factorization interpolation, done in routine TRID*FS
! ###################################################################
! -------------------------------------------------------------------
! Periodic case; pentadiagonal
! -------------------------------------------------------------------
    if ((stagger_on) .and. g%periodic) then
        g%lu0i => x(:, ig:)

        select case (g%mode_fdm)
        case DEFAULT
            call FDM_C0INT6P_LHS(nx, g%lu0i(1, 1), g%lu0i(1, 2), g%lu0i(1, 3))
        end select
        call TRIDPFS(nx, g%lu0i(1, 1), g%lu0i(1, 2), g%lu0i(1, 3), g%lu0i(1, 4), g%lu0i(1, 5))
        ig = ig + 5
    end if

! ###################################################################
! LU factorization first interp. derivative, done in routine TRID*FS
! ###################################################################
! -------------------------------------------------------------------
! Periodic case; pentadiagonal
! -------------------------------------------------------------------
    if ((stagger_on) .and. g%periodic) then
        g%lu1i => x(:, ig:)

        select case (g%mode_fdm)
        case DEFAULT
            call FDM_C1INT6P_LHS(nx, g%jac, g%lu1i(1, 1), g%lu1i(1, 2), g%lu1i(1, 3))
        end select
        call TRIDPFS(nx, g%lu1i(1, 1), g%lu1i(1, 2), g%lu1i(1, 3), g%lu1i(1, 4), g%lu1i(1, 5))
        ig = ig + 5
    end if

! ###################################################################
! Modified wavenumbers in periodic case
! ###################################################################
    g%mwn => x(:, ig:)

    if (g%periodic) then
        do i = 1, nx ! Define wavenumbers
            if (i <= nx/2 + 1) then
                wrk1d(i, 1) = 2.0_wp*pi_wp*real(i - 1, wp)/real(nx, wp)
            else
                wrk1d(i, 1) = 2.0_wp*pi_wp*real(i - 1 - nx, wp)/real(nx, wp)
            end if
        end do

! -------------------------------------------------------------------
! First-order derivative
! -------------------------------------------------------------------
        if (.not. stagger_on) then

            select case (g%mode_fdm)

            case (FDM_COM6_JACOBIAN, FDM_COM6_DIRECT)
                a1 = 1.0_wp/3.0_wp
                a2 = 0.0_wp
                b1 = 7.0_wp/9.0_wp
                b2 = 1.0_wp/36.0_wp
                b3 = 0.0_wp

            case (FDM_COM6_JACPENTA)
                a1 = C1N6M_ALPHA2/2.0_wp
                a2 = C1N6M_BETA2/2.0_wp
                b1 = C1N6M_A/2.0_wp
                b2 = C1N6M_BD2/2.0_wp
                b3 = C1N6M_CD3/2.0_wp

            case (FDM_COM8_JACOBIAN)
                a1 = 3.0_wp/8.0_wp
                a2 = 0.0_wp
                b1 = 25.0_wp/32.0_wp
                b2 = 1.0_wp/20.0_wp
                b3 = -1.0_wp/480.0_wp

            end select

            g%mwn(:, 1) = 2.0_wp*(b1*sin(wrk1d(:, 1)) + b2*sin(2.0_wp*wrk1d(:, 1)) + b3*sin(3.0_wp*wrk1d(:, 1))) &
                          /(1.0_wp + 2.0_wp*a1*cos(wrk1d(:, 1)) + 2.0_wp*a2*cos(wrk1d(:, 1)))

        else ! staggered case has different modified wavenumbers!

            select case (g%mode_fdm)

            case DEFAULT
                a1 = 9.0_wp/62.0_wp
                b1 = 63.0_wp/62.0_wp
                b2 = 17.0_wp/62.0_wp

            end select

            g%mwn(:, 1) = 2.0_wp*(b1*sin(1.0_wp/2.0_wp*wrk1d(:, 1)) + b2/3.0_wp*sin(3.0_wp/2.0_wp*wrk1d(:, 1))) &
                          /(1.0_wp + 2.0_wp*a1*cos(wrk1d(:, 1)))

        end if

! Final calculations because it is mainly used in the Poisson solver like this
        g%mwn(:, 1) = (g%mwn(:, 1)/g%jac(1, 1))**2

! -------------------------------------------------------------------
! Second-order derivative
! -------------------------------------------------------------------
        select case (g%mode_fdm)

        case (FDM_COM4_JACOBIAN) ! Not yet implemented

        ! case (FDM_COM6_DIRECT, FDM_COM6_JACPENTA)
        !     a1 = 2.0_wp/11.0_wp         ! Lele's standard 6th-order pentadiagonal compact
        !     b1 = 12.0_wp/11.0_wp
        !     b2 = 3.0_wp/44.0_wp
        !     b3 = 0.0_wp
            
        case (FDM_COM6_JACOBIAN, FDM_COM6_DIRECT, FDM_COM6_JACPENTA)
            kc = pi_wp**2.0_wp          ! Lambellais' 6th-order heptadiagonal compact
            a1 = (272.0_wp - 45.0_wp*kc)/(416.0_wp - 90.0_wp*kc)
            b1 = (48.0_wp - 135.0_wp*kc)/(1664.0_wp - 360.0_wp*kc)
            b2 = (528.0_wp - 81.0_wp*kc)/(208.0_wp - 45.0_wp*kc)/4.0_wp
            b3 = -(432.0_wp - 63.0_wp*kc)/(1664.0_wp - 360.0_wp*kc)/9.0_wp

        case (FDM_COM8_JACOBIAN) ! Not yet implemented

        end select

        g%mwn(:, 2) = 2.0_wp*(b1*(1.0_wp - cos(wrk1d(:, 1))) + b2*(1.0_wp - cos(2.0_wp*wrk1d(:, 1))) + b3*(1.0_wp - cos(3.0_wp*wrk1d(:, 1)))) &
                      /(1.0_wp + 2.0_wp*a1*cos(wrk1d(:, 1)))

! Final calculations because it is mainly used in the Helmholtz solver like this
        g%mwn(:, 2) = g%mwn(:, 2)/(g%jac(1, 1)**2)

        ig = ig + 2

    end if

! ###################################################################
! Density correction in anelastic mode
! ###################################################################
    g%rhoinv => x(:, ig)

    g%anelastic = .false. ! Default; activated in FI_BACKGROUND_INITIALIZE

    ig = ig + 1
    
! ###################################################################
! Check array sizes
! ###################################################################
    ! IF ( ig .NE. g%inb_grid ) THEN
    !    CALL TLAB_WRITE_ASCII(efile, 'FDM_INITIALIZE. Grid size incorrect.')
    !    CALL TLAB_STOP(DNS_ERROR_DIMGRID)
    ! ENDIF

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'Leaving SUBOURINTE FDM_INITIALIZE')
#endif

    return
end subroutine FDM_INITIALIZE
