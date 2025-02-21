
program VPARTIAL2
    use TLab_Constants, only: wp, wi
    use FDM, only: fdm_dt
    use OPR_PARTIAL
    implicit none

    type(fdm_dt) :: g
    integer(wi) imax, jmax, kmax, i, wk, idummy
    parameter(imax=128)
    real(wp) x(imax, 3 + 4*3 + 4*3), u(imax), du1(imax), du2(imax), due(imax)
    real(wp) wrk1d(imax, 5), wrk2d(imax), wrk3d(imax), bcs(2, 2)
    real(wp) tmp(imax)

! ###################################################################
    g%size = imax
    g%scale = C_2_R*C_PI_R
    g%mode_fdm1 = FDM_COM6_JACOBIAN
    g%uniform = .true.
    jmax = 1
    kmax = 1

    write (*, *) 'Periodic (.TRUE. or .FALSE.)?'
    read (*, *) g%periodic
    write (*, *) 'Wavenumber ?'
    read (*, *) wk

! CHANGE TO UPDATE NEW fdm_dt
    if (g%periodic) then
        do i = 1, imax
            g%nodes(i) = M_REAL(i - 1)/M_REAL(imax)*g%scale
        end do
    else
        open (21, file='y.dat')
        do i = 1, imax
!        g%nodes(i) = M_REAL(i-1)/M_REAL(imax-1)*g%scale
            read (21, *) idummy, g%nodes(i)
        end do
        close (21)
        g%scale = g%nodes(imax) - g%nodes(1)
    end if

    call FDM_Initialize(g, wrk1d)

! ###################################################################
! Define the function
    do i = 1, imax
        u(i) = sin(C_2_R*C_PI_R/g%scale*M_REAL(wk)*g%nodes(i))
        due(i) = -(C_2_R*C_PI_R/g%scale*M_REAL(wk))**2*u(i)
!     u(i) = EXP(-(g%nodes(i)-C_PI_R)**2/(C_PI_R**2/64.))
!     due(i) = -C_2_R*(g%nodes(i)-C_PI_R)/(C_PI_R**2/64.)*u(i)
    end do

! ###################################################################
    bcs(:, 1) = 0
    bcs(:, 2) = 1
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs(1, 1), g, u, tmp)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs(1, 2), g, tmp, du1)
    call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs, g, u, du2, tmp)

    open (20, file='partial.dat')
    do i = 1, imax
        write (20, '(7e)') g%nodes(i), g%jac(i, 1), g%jac(i, 2), u(i), due(i), du1(i), du2(i)
    end do
    close (20)

    stop
end program VPARTIAL2
