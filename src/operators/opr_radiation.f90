#include "types.h"
#include "dns_const.h"

!########################################################################
!#
!# Calculate the radiation term assuming the 1D approximation and
!# using compact schemes to calculate the integral term.
!#
!########################################################################
subroutine OPR_RADIATION(radiation, nx, ny, nz, g, s, r)
    use TLAB_TYPES, only: term_dt, grid_dt
    use TLAB_ARRAYS, only: wrk1d, wrk3d
    use AVGS, only: AVG1V2D_V

    implicit none

    type(term_dt), intent(IN) :: radiation
    TINTEGER, intent(IN) :: nx, ny, nz
    type(grid_dt), intent(IN) :: g
    TREAL, dimension(nx*ny*nz), intent(IN) :: s        ! Radiatively active scalar
    TREAL, dimension(nx*ny*nz), intent(OUT) :: r        ! Radiative heating rate

    target s, r

! -----------------------------------------------------------------------
    TINTEGER j, ip, ip2, nxy, nxz, ibc
    TREAL delta_inv, f0, f1
    TREAL, dimension(:), pointer :: p_org, p_dst

! #######################################################################
    nxy = nx*ny ! For transposition to make y direction the last one
    ibc = 2     ! Boundary condition at the top for integral calulation

! Prepare the pentadiagonal system
    call INT_C1N6_LHS(ny, ibc, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
    call PENTADFS(ny - 1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))

    delta_inv = C_1_R/radiation%parameters(2)

! ###################################################################
    if (radiation%type == EQNS_RAD_BULK1D_GLOBAL) then
        nxz = 1

        if (nz == 1) then
            p_org => wrk3d
            p_dst => r
        else
            p_org => r
            p_dst => wrk3d
        end if

        call AVG1V2D_V(nx, ny, nz, 1, s, p_org, p_dst) ! Calculate averaged scalar into p_org; p_dst is auxiliar

! ###################################################################
    else if (radiation%type == EQNS_RAD_BULK1D_LOCAL) then
        nxz = nx*nz

        if (nz == 1) then
            p_org => s
            p_dst => r
        else
            p_org => r
            p_dst => wrk3d

#ifdef USE_ESSL
            call DGETMO(s, nxy, nxy, nz, p_org, nz)
#else
            call DNS_TRANSPOSE(s, nxy, nz, nxy, p_org, nz)
#endif
        end if

    end if

! ###################################################################
! Calculate (negative) integral path.
    call INT_C1N6_RHS(ny, nxz, ibc, g%jac, p_org, p_dst)
    call PENTADSS(ny - 1, nxz, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), p_dst)
    ip = nxz*(ny - 1) + 1; p_dst(ip:ip + nxz - 1) = C_0_R ! boundary condition

! Calculate radiative heating rate
    f0 = radiation%parameters(1)
    if (ABS(radiation%parameters(3)) > C_0_R) then
        f1 = radiation%parameters(3)
        do j = ny, 1, -1
            ip = nx*nz*(j - 1) + 1
            ip2 = ip + nx*nz - 1
            p_dst(ip:ip2) = p_org(ip:ip2)*( &
                            EXP(p_dst(ip:ip2)*delta_inv)*f0 &
                            + EXP((p_dst(1:nx*nz) - p_dst(ip:ip2))*delta_inv)*f1)
        end do
    else
        do j = 1, ny*nxz
            p_dst(j) = p_org(j)*EXP(p_dst(j)*delta_inv)*f0
        end do
!  p_dst(1:ny*nxz) = radiation%parameters(1) *p_org(1:ny*nxz) *DEXP( p_dst(1:ny*nxz) *delta_inv ) seg-fault; need ulimit -u unlimited
    end if

! ###################################################################
    if (radiation%type == EQNS_RAD_BULK1D_GLOBAL) then
        do j = ny, 1, -1
            ip = nx*nz*(j - 1) + 1; p_dst(ip:ip + nx*nz - 1) = p_dst(j)
        end do
    end if

    if (nz > 1) then ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(p_dst, nz, nz, nxy, r, nxy)
#else
        call DNS_TRANSPOSE(p_dst, nz, nxy, nz, r, nxy)
#endif
    end if

    nullify (p_org, p_dst)

    return
end subroutine OPR_RADIATION

!########################################################################
!########################################################################
subroutine OPR_RADIATION_FLUX(radiation, nx, ny, nz, g, s, r)
    use TLAB_TYPES, only: term_dt, grid_dt
    use TLAB_ARRAYS, only: wrk1d, wrk3d
    use AVGS, only: AVG1V2D_V

    implicit none

    type(term_dt), intent(IN) :: radiation
    TINTEGER, intent(IN) :: nx, ny, nz
    type(grid_dt), intent(IN) :: g
    TREAL, dimension(nx*ny*nz), intent(IN) :: s        ! Radiatively active scalar
    TREAL, dimension(nx*ny*nz), intent(OUT) :: r        ! Radiative flux

    target s, r

! -----------------------------------------------------------------------
    TINTEGER j, ip, ip2, nxy, nxz, ibc
    TREAL delta_inv, f0, f1
    TREAL, dimension(:), pointer :: p_org, p_dst

! #######################################################################
    nxy = nx*ny ! For transposition to make y direction the last one
    ibc = 2     ! Boundary condition at the top for integral calulation

! Prepare the pentadiagonal system
    call INT_C1N6_LHS(ny, ibc, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
    call PENTADFS(ny - 1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))

    delta_inv = C_1_R/radiation%parameters(2)

! ###################################################################
    if (radiation%type == EQNS_RAD_BULK1D_GLOBAL) then
        nxz = 1

        if (nz == 1) then
            p_org => wrk3d
            p_dst => r
        else
            p_org => r
            p_dst => wrk3d
        end if

        call AVG1V2D_V(nx, ny, nz, 1, s, p_org, p_dst) ! Calculate averaged scalar into p_org; p_dst is auxiliar

! ###################################################################
    else if (radiation%type == EQNS_RAD_BULK1D_LOCAL) then
        nxz = nx*nz

        if (nz == 1) then
            p_org => s
            p_dst => r
        else
            p_org => r
            p_dst => wrk3d

#ifdef USE_ESSL
            call DGETMO(s, nxy, nxy, nz, p_org, nz)
#else
            call DNS_TRANSPOSE(s, nxy, nz, nxy, p_org, nz)
#endif
        end if

    end if

! ###################################################################
! Calculate (negative) integral path.
    call INT_C1N6_RHS(ny, nxz, ibc, g%jac, p_org, p_dst)
    call PENTADSS(ny - 1, nxz, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), p_dst)
    ip = nxz*(ny - 1) + 1; p_dst(ip:ip + nxz - 1) = C_0_R ! boundary condition

! Calculate radiative flux
    f0 = -radiation%parameters(1)*radiation%parameters(2)
    if (ABS(radiation%parameters(3)) > C_0_R) then
        f1 = -radiation%parameters(3)*radiation%parameters(2)
        do j = ny, 1, -1
            ip = nx*nz*(j - 1) + 1
            ip2 = ip + nx*nz - 1
            p_dst(ip:ip2) = EXP(p_dst(ip:ip2)*delta_inv)*f0 &
                            - EXP((p_dst(1:nx*nz) - p_dst(ip:ip2))*delta_inv)*f1
        end do
    else
        do j = 1, ny*nxz
            p_dst(j) = EXP(p_dst(j)*delta_inv)*f0
        end do
!  p_dst(1:ny*nxz) = radiation%parameters(1) *p_org(1:ny*nxz) *DEXP( p_dst(1:ny*nxz) *delta_inv ) seg-fault; need ulimit -u unlimited
    end if

! ###################################################################
    if (radiation%type == EQNS_RAD_BULK1D_GLOBAL) then
        do j = ny, 1, -1
            ip = nx*nz*(j - 1) + 1; p_dst(ip:ip + nx*nz - 1) = p_dst(j)
        end do
    end if

    if (nz > 1) then ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(p_dst, nz, nz, nxy, r, nxy)
#else
        call DNS_TRANSPOSE(p_dst, nz, nxy, nz, r, nxy)
#endif
    end if

    nullify (p_org, p_dst)

    return
end subroutine OPR_RADIATION_FLUX
