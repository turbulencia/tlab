#include "dns_const.h"
#include "dns_error.h"

module RADIATION_M
    use TLAB_CONSTANTS, only: wp, wi, BCS_MAX, BCS_MIN, efile, MAX_PROF
    use TLAB_TYPES, only: term_dt, grid_dt
    use TLAB_VARS, only: imode_eqns, inb_scal_array
    use TLAB_VARS, only: radiation
    use TLAB_ARRAYS, only: wrk3d
    use TLAB_PROCS, only: TLAB_WRITE_ASCII, TLAB_STOP
    use THERMO_VARS, only: imixture
    use OPR_PARTIAL, only: OPR_PARTIAL_Y
    use OPR_ODES
    implicit none

    integer, parameter :: TYPE_NONE = 0
    integer, parameter :: TYPE_LW_BULK1D_LIQUID = 1
    integer, parameter :: TYPE_LW_BULK1D = 2

    real(wp) :: sigma = 5.67037442e-8_wp ! W /m^2 /K

    public :: RADIATION_INITIALIZE
    public :: RADIATION_X

contains
    ! ###################################################################
    ! ###################################################################
    subroutine RADIATION_READBLOCK(bakfile, inifile, block, var)
        character(len=*), intent(in) :: bakfile, inifile, block
        type(term_dt), intent(out) :: var

        character(len=512) sRes
        integer(wi) idummy

        ! -------------------------------------------------------------------
        call TLAB_WRITE_ASCII(bakfile, '#')
        call TLAB_WRITE_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLAB_WRITE_ASCII(bakfile, '#Type=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Scalar=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Scalar=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

        if (var%type == EQNS_NONE) then        ! backwards compatibility, to be removed
            call SCANINICHAR(bakfile, inifile, block, 'Type', 'None', sRes)
            if (trim(adjustl(sRes)) == 'none') then; var%type = TYPE_NONE
            else if (trim(adjustl(sRes)) == 'irbulk1dliquid') then; var%type = TYPE_LW_BULK1D_LIQUID
            else if (trim(adjustl(sRes)) == 'irbulk1d') then; var%type = TYPE_LW_BULK1D
            else if (trim(adjustl(sRes)) == 'bulk1dlocal') then; var%type = EQNS_RAD_BULK1D_LOCAL
            else if (trim(adjustl(sRes)) == 'bulk1dglobal') then; var%type = EQNS_RAD_BULK1D_GLOBAL
            else
                call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong Radiation option.')
                call TLAB_STOP(DNS_ERROR_OPTION)
            end if
        end if

        var%active = .false.
        if (var%type /= EQNS_NONE) then
            call SCANINIINT(bakfile, inifile, block, 'Scalar', '1', idummy)
            var%active(idummy) = .true.

            var%parameters(:) = 0.0_wp
            call SCANINICHAR(bakfile, inifile, block, 'Parameters', '1.0', sRes)
            idummy = MAX_PROF
            call LIST_REAL(sRes, idummy, var%parameters)

        end if

        return
    end subroutine RADIATION_READBLOCK

!########################################################################
!########################################################################
    subroutine RADIATION_INITIALIZE()

        ! -------------------------------------------------------------------
        ! By default, transport and radiation are caused by last scalar
        radiation%scalar = inb_scal_array

        if (imixture == MIXT_TYPE_AIRWATER .or. imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            if (radiation%type /= EQNS_NONE) then
                radiation%active(inb_scal_array) = .true. ! liquid
                radiation%active(inb_scal_array + 1) = .true. ! buoyancy
            end if

        end if

        ! -------------------------------------------------------------------
        ! in case nondimensional we need to adjust sigma

        return
    end subroutine RADIATION_INITIALIZE

!########################################################################
!########################################################################
    subroutine IR_BULK1_LIQUID(radiation, nx, ny, nz, g, a, source, flux)
        type(term_dt), intent(in) :: radiation
        integer(wi), intent(in) :: nx, ny, nz
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: a(nx*nz, ny)                ! bulk absorption coefficent
        real(wp), intent(out) :: source(nx*nz, ny)
        real(wp), intent(out), optional :: flux(nx*nz, ny)

        target a, source, flux

! -----------------------------------------------------------------------
        integer(wi) j, nxy, nxz
        real(wp) f0, f1
        real(wp), pointer :: p_org(:, :) => null()
        real(wp), pointer :: p_tau(:, :) => null()
        real(wp), pointer :: p_source(:, :) => null()
        real(wp), pointer :: p_flux(:, :) => null()

! #######################################################################
        nxy = nx*ny     ! For transposition to make y direction the last one
        nxz = nx*nz

        if (nz == 1) then
            p_org => a
            p_tau(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
            p_source => source
            if (present(flux)) p_flux => flux
        else
            p_org => source
            p_tau(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
            if (present(flux)) then
                p_source => flux
                p_flux(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
            else
                p_source(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
            end if

#ifdef USE_ESSL
            call DGETMO(a, nxy, nxy, nz, p_org, nz)
#else
            call DNS_TRANSPOSE(a, nxy, nz, nxy, p_org, nz)
#endif
        end if

! ###################################################################
        ! Calculate (negative) optical thickness; integrating from the top, to be checked
        p_tau(:, ny) = 0.0_wp   ! boundary condition
        call OPR_Integral1(nxz, g, p_org, p_tau, BCS_MAX)
        ! p_tau(:, 1) = 0.0_wp     ! boundary condition
        ! call OPR_Integral1(nxz, g, p_org, p_tau, BCS_MIN)

! Calculate exp(-\tau(z, zmax)/\mu)
        do j = ny, 1, -1
            p_tau(:, j) = exp(p_tau(:, j))
        end do
        !  p_tau = dexp(p_tau)         seg-fault; need ulimit -u unlimited

! ###################################################################
! Calculate heating rate
        f0 = radiation%parameters(1)*radiation%parameters(2)
        f1 = radiation%parameters(3)*radiation%parameters(2)
        if (abs(radiation%parameters(3)) > 0.0_wp) then
            do j = ny, 1, -1
                p_source(:, j) = p_org(:, j)*(p_tau(:, j)*f0 &                       ! downward flux
                                              + p_tau(:, 1)/p_tau(:, j)*f1)       ! upward flux
            end do
        else
            p_source = p_org*p_tau*f0
        end if

        if (nz > 1) then ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
            call DGETMO(p_source, nz, nz, nxy, source, nxy)
#else
            call DNS_TRANSPOSE(p_source, nz, nxy, nz, source, nxy)
#endif
        end if

! ###################################################################
! Calculate radiative flux, if necessary
        if (present(flux)) then
            f0 = -radiation%parameters(1)*radiation%parameters(2)
            f1 = radiation%parameters(3)*radiation%parameters(2)
            if (abs(radiation%parameters(3)) > 0.0_wp) then
                do j = ny, 1, -1
                    p_flux(:, j) = p_tau(:, j)*f0 &                       ! downward flux
                                   + p_tau(:, 1)/p_tau(:, j)*f1       ! upward flux
                end do
            else
                p_flux = p_tau*f0
            end if

            if (nz > 1) then ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
                call DGETMO(p_flux, nz, nz, nxy, flux, nxy)
#else
                call DNS_TRANSPOSE(p_flux, nz, nxy, nz, flux, nxy)
#endif
            end if

        end if

! ###################################################################
        nullify (p_org, p_tau, p_source, p_flux)

        return
    end subroutine IR_BULK1_LIQUID

!########################################################################
!########################################################################
    subroutine RADIATION_X(radiation, nx, ny, nz, g, s, source, a, b, flux)
        use THERMO_ANELASTIC
        type(term_dt), intent(in) :: radiation
        integer(wi), intent(in) :: nx, ny, nz
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: s(nx*ny*nz, inb_scal_array)
        real(wp), intent(out) :: source(nx*ny*nz)
        real(wp), intent(inout) :: a(nx*ny*nz)      ! space for bulk absorption coefficient
        real(wp), intent(inout) :: b(nx*ny*nz)      ! space for emission function
        real(wp), intent(out), optional :: flux(nx*ny*nz)

        ! -----------------------------------------------------------------------
        real(wp) delta_inv

        !########################################################################
        select case (radiation%type)
        case (TYPE_LW_BULK1D_LIQUID, EQNS_RAD_BULK1D_LOCAL)
            ! bulk absorption coefficient
            delta_inv = 1.0_wp/radiation%parameters(2)
            if (imode_eqns == DNS_EQNS_ANELASTIC) then
                call THERMO_ANELASTIC_WEIGHT_OUTPLACE(nx, ny, nz, rbackground, s(:, radiation%scalar(1)), a)
                a = delta_inv*a
            else
                a = delta_inv*s(:, radiation%scalar(1))
            end if

            if (present(flux)) then
                call IR_BULK1_LIQUID(radiation, nx, ny, nz, g, a, source, flux)
            else
                call IR_BULK1_LIQUID(radiation, nx, ny, nz, g, a, source)
            end if

        case (TYPE_LW_BULK1D)
            ! bulk absorption coefficient

            ! emission function

            ! if (present(flux)) then
            !     call IR_BULK1(radiation, nx, ny, nz, g, a, b, source, flux)
            ! else
            !     call IR_BULK1(radiation, nx, ny, nz, g, a, b, source)
            ! end if

        end select

    end subroutine RADIATION_X

end module
