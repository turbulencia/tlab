#include "dns_const.h"

!########################################################################
!#
!# Calculate the radiation term assuming the 1D approximation and
!# using compact schemes to calculate the integral term.
!#
!########################################################################
subroutine OPR_RADIATION(radiation, nx, ny, nz, g, s, r)
    use TLAB_CONSTANTS, only: wp, wi, BCS_MAX!, BCS_MIN
    use TLAB_TYPES, only: term_dt, grid_dt
    use TLAB_ARRAYS, only: wrk3d
    use AVGS, only: AVG1V2D_V
    use OPR_ODES

    implicit none

    type(term_dt), intent(IN) :: radiation
    integer(wi), intent(IN) :: nx, ny, nz
    type(grid_dt), intent(IN) :: g
    real(wp), intent(IN) :: s(nx*ny*nz)         ! Radiatively active scalar
    real(wp), intent(OUT) :: r(nx*ny*nz)        ! Radiative heating rate

    target s, r

! -----------------------------------------------------------------------
    integer(wi) j, ip, ip2, nxy, nxz
    real(wp) delta_inv, f0, f1
    real(wp), pointer :: p_org(:), p_dst(:)

! #######################################################################
    nxy = nx*ny ! For transposition to make y direction the last one

    delta_inv = 1.0_wp/radiation%parameters(2)

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
! Calculate (negative) integral path; integrating from the top, to be checked
    ip = nxz*(ny - 1) + 1; p_dst(ip:ip + nxz - 1) = 0.0_wp ! boundary condition
    call OPR_Integral1(nxz, g, p_org, p_dst, BCS_MAX)
    ! p_dst(1:nxz) = 0.0_wp ! boundary condition
    ! call OPR_Integral1(nxz, g, p_org, p_dst, BCS_MIN)

! Calculate radiative heating rate
    f0 = radiation%parameters(1)
    if (ABS(radiation%parameters(3)) > 0.0_wp) then
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
    use TLAB_CONSTANTS, only: wp, wi, BCS_MAX!, BCS_MIN
    use TLAB_TYPES, only: term_dt, grid_dt
    use TLAB_ARRAYS, only: wrk3d
    use AVGS, only: AVG1V2D_V
    use OPR_ODES

    implicit none

    type(term_dt), intent(IN) :: radiation
    integer(wi), intent(IN) :: nx, ny, nz
    type(grid_dt), intent(IN) :: g
    real(wp), dimension(nx*ny*nz), intent(IN) :: s          ! Radiatively active scalar
    real(wp), dimension(nx*ny*nz), intent(OUT) :: r         ! Radiative flux

    target s, r

! -----------------------------------------------------------------------
    integer(wi) j, ip, ip2, nxy, nxz
    real(wp) delta_inv, f0, f1
    real(wp), dimension(:), pointer :: p_org, p_dst

! #######################################################################
    nxy = nx*ny ! For transposition to make y direction the last one

    delta_inv = 1.0_wp/radiation%parameters(2)

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
! Calculate (negative) integral path; integrating from the top, to be checked
    ip = nxz*(ny - 1) + 1; p_dst(ip:ip + nxz - 1) = 0.0_wp ! boundary condition
    call OPR_Integral1(nxz, g, p_org, p_dst, BCS_MAX)
    ! p_dst(1:nxz) = 0.0_wp ! boundary condition
    ! call OPR_Integral1(nxz, g, p_org, p_dst, BCS_MIN)

! Calculate radiative flux
    f0 = -radiation%parameters(1)*radiation%parameters(2)
    if (ABS(radiation%parameters(3)) > 0.0_wp) then
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
