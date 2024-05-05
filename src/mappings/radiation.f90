#include "dns_const.h"
#include "dns_error.h"

module Radiation
    use TLAB_CONSTANTS, only: wp, wi, pi_wp, BCS_MAX, BCS_MIN, efile, MAX_PROF
    use TLAB_TYPES, only: term_dt, grid_dt
    use TLAB_VARS, only: imode_eqns, inb_scal_array
    ! use TLAB_VARS, only: infrared
    use TLAB_ARRAYS, only: wrk3d
    use TLAB_PROCS, only: TLAB_WRITE_ASCII, TLAB_STOP
    use THERMO_VARS, only: imixture
    use OPR_PARTIAL, only: OPR_PARTIAL_Y
    use OPR_ODES
    implicit none

    integer, parameter :: TYPE_NONE = 0
    integer, parameter :: TYPE_IR_BULK1D_LIQUID = 1
    integer, parameter :: TYPE_IR_BULK1D = 2
    integer, parameter :: TYPE_BULK1DLOCAL = 10         ! backwards compatibility, to be removed

    real(wp), parameter :: sigma = 5.67037442e-8_wp     ! W /m^2 /K
    real(wp) :: sigma_o_pi

    public :: Radiation_Initialize
    public :: Radiation_Infrared

contains
!########################################################################
!########################################################################
    subroutine Radiation_Initialize(inifile)
        use TLAB_VARS, only: infrared
        character(len=*), intent(in) :: inifile

        character(len=32) bakfile
        character(len=512) sRes
        integer(wi) idummy

        ! -------------------------------------------------------------------
        bakfile = trim(adjustl(inifile))//'.bak'

        call TLAB_WRITE_ASCII(bakfile, '#')
        call TLAB_WRITE_ASCII(bakfile, '#[Radiation]')
        call TLAB_WRITE_ASCII(bakfile, '#Type=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Scalar=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

        call SCANINICHAR(bakfile, inifile, 'Radiation', 'Type', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') &
            call SCANINICHAR(bakfile, inifile, 'Main', 'TermRadiation', 'None', sRes)               ! backwards compatibility, to be removed
        if (trim(adjustl(sRes)) == 'none') then; infrared%type = TYPE_NONE
        else if (trim(adjustl(sRes)) == 'irbulk1dliquid') then; infrared%type = TYPE_IR_BULK1D_LIQUID
        else if (trim(adjustl(sRes)) == 'irbulk1d') then; infrared%type = TYPE_IR_BULK1D
        else if (trim(adjustl(sRes)) == 'bulk1dlocal') then; infrared%type = TYPE_BULK1DLOCAL    ! backwards compatibility, to be removed
        else
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong Radiation option.')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if

        infrared%active = .false.
        if (infrared%type /= TYPE_NONE) then
            call SCANINIINT(bakfile, inifile, 'Radiation', 'Scalar', '1', idummy)
            infrared%active(idummy) = .true.

            infrared%parameters(:) = 0.0_wp
            call SCANINICHAR(bakfile, inifile, 'Radiation', 'Parameters', '1.0', sRes)
            idummy = MAX_PROF
            call LIST_REAL(sRes, idummy, infrared%parameters)

        end if

        ! -------------------------------------------------------------------
        ! By default, transport and radiation are caused by last scalar
        infrared%scalar = inb_scal_array

        if (imixture == MIXT_TYPE_AIRWATER .or. imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            if (infrared%type /= EQNS_NONE) then
                infrared%active(inb_scal_array) = .true. ! liquid
                infrared%active(inb_scal_array + 1) = .true. ! buoyancy
            end if

        end if

        ! backwards compatibility
        if (infrared%type == TYPE_BULK1DLOCAL) then
            infrared%parameters(1) = infrared%parameters(1)*infrared%parameters(2)
            infrared%parameters(3) = infrared%parameters(3)*infrared%parameters(2)
            infrared%parameters(2) = 1.0_wp/infrared%parameters(2)
        end if

        ! -------------------------------------------------------------------
        ! in case nondimensional we need to adjust sigma

        sigma_o_pi = sigma/pi_wp
        ! sigma_o_pi = 0.0_wp       ! testing

        return
    end subroutine Radiation_Initialize

!########################################################################
!########################################################################
    subroutine IR_Bulk1D_Liquid(infrared, nx, ny, nz, g, a_source, flux)
        type(term_dt), intent(in) :: infrared
        integer(wi), intent(in) :: nx, ny, nz
        type(grid_dt), intent(in) :: g
        real(wp), intent(inout) :: a_source(nx*nz, ny)      ! input as bulk absorption coefficent, output as source
        real(wp), intent(out), optional :: flux(nx*nz, ny)

        target a_source, flux

! -----------------------------------------------------------------------
        integer(wi) j, nxy, nxz
        real(wp) fd, fu
        real(wp), pointer :: p_a(:, :) => null()
        real(wp), pointer :: p_tau(:, :) => null()
        real(wp), pointer :: p_source(:, :) => null()
        real(wp), pointer :: p_flux(:, :) => null()

! #######################################################################
        nxy = nx*ny     ! For transposition to make y direction the last one
        nxz = nx*nz

        if (nz == 1) then
            p_a => a_source
            p_source => a_source
            p_tau(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
            if (present(flux)) p_flux => flux
        else
            p_a(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
            p_source(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
            if (present(flux)) then
                p_tau => flux
                p_flux(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
            else
                p_tau => a_source
            end if

#ifdef USE_ESSL
            call DGETMO(a_source, nxy, nxy, nz, p_a, nz)
#else
            call DNS_TRANSPOSE(a_source, nxy, nz, nxy, p_a, nz)
#endif
        end if

! ###################################################################
        ! Calculate (negative) optical thickness; integrating from the top, to be checked
        p_tau(:, ny) = 0.0_wp   ! boundary condition
        call OPR_Integral1(nxz, g, p_a, p_tau, BCS_MAX)
        ! p_tau(:, 1) = 0.0_wp     ! boundary condition
        ! call OPR_Integral1(nxz, g, p_org, p_tau, BCS_MIN)

! Calculate exp(-tau(z, zmax)/\mu)
        do j = ny, 1, -1
            p_tau(:, j) = exp(p_tau(:, j))
        end do
        !  p_tau = dexp(p_tau)         seg-fault; need ulimit -u unlimited

! ###################################################################
! Calculate heating rate
        fd = infrared%parameters(1)
        fu = infrared%parameters(3)
        if (abs(infrared%parameters(3)) > 0.0_wp) then
            do j = ny, 1, -1
                p_source(:, j) = p_a(:, j)*(p_tau(:, j)*fd &                       ! downward flux
                                            + p_tau(:, 1)/p_tau(:, j)*fu)       ! upward flux
            end do
        else
            p_source = p_a*p_tau*fd
        end if

        if (nz > 1) then ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
            call DGETMO(p_source, nz, nz, nxy, a_source, nxy)
#else
            call DNS_TRANSPOSE(p_source, nz, nxy, nz, a_source, nxy)
#endif
        end if

! ###################################################################
! Calculate radiative flux, if necessary
        if (present(flux)) then
            fd = -infrared%parameters(1)
            fu = infrared%parameters(3)
            if (abs(infrared%parameters(3)) > 0.0_wp) then
                do j = ny, 1, -1
                    p_flux(:, j) = fd*p_tau(:, j) &                       ! downward flux
                                   + fu*p_tau(:, 1)/p_tau(:, j)       ! upward flux
                end do
            else
                p_flux = p_tau*fd
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
        nullify (p_a, p_tau, p_source, p_flux)

        return
    end subroutine IR_Bulk1D_Liquid

!########################################################################
!########################################################################
! Here we do not treat separately 2d and 3d cases for the trasposition because it was a bit messy...
    subroutine IR_Bulk1D(infrared, nx, ny, nz, g, a_source, b, tmp1, tmp2, flux)
        type(term_dt), intent(in) :: infrared
        integer(wi), intent(in) :: nx, ny, nz
        type(grid_dt), intent(in) :: g
        real(wp), intent(inout) :: a_source(nx*nz, ny)          ! input as bulk absorption coefficent, output as source
        real(wp), intent(inout) :: b(nx*nz, ny)                 ! emission function
        real(wp), intent(inout) :: tmp1(nx*nz, ny), tmp2(nx*nz, ny)
        real(wp), intent(out), optional :: flux(nx*nz, ny)

        target a_source, b, tmp1

! -----------------------------------------------------------------------
        integer(wi) j, nxy, nxz
        real(wp) fd, mu, epsilon, dummy!, fu
        real(wp), pointer :: p_a(:, :) => null()
        real(wp), pointer :: p_tau(:, :) => null()
        real(wp), pointer :: p_ab(:, :) => null()
        real(wp), pointer :: p_source(:, :) => null()
        real(wp), pointer :: p_flux(:, :) => null()

        real(wp) :: bcs(nx*nz)

! #######################################################################
        nxy = nx*ny     ! For transposition to make y direction the last one
        nxz = nx*nz

        p_a(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
        p_source(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
        p_ab => a_source
        p_tau => b
        p_flux => tmp1

#ifdef USE_ESSL
        call DGETMO(a_source, nxy, nxy, nz, p_a, nz)
        call DGETMO(b, nxy, nxy, nz, p_ab, nz)
#else
        call DNS_TRANSPOSE(a_source, nxy, nz, nxy, p_a, nz)
        call DNS_TRANSPOSE(b, nxy, nz, nxy, p_ab, nz)
#endif

! ###################################################################
        mu = 0.5_wp*(1.0_wp/sqrt(3.0_wp) + 1.0_wp/sqrt(2.0_wp))     ! solid angle factor, in (1/sqrt{3}, 1/sqrt{2})
        ! mu = 1.0_wp     ! testing

        dummy = 1.0_wp/mu
        p_a = p_a*dummy
        ! Calculate (negative) optical thickness; integrating from the top, to be checked
        p_tau(:, ny) = 0.0_wp   ! boundary condition
        call OPR_Integral1(nxz, g, p_a, p_tau, BCS_MAX)
        ! p_tau(:, 1) = 0.0_wp     ! boundary condition
        ! call OPR_Integral1(nxz, g, p_org, p_tau, BCS_MIN)

        ! Calculate exp(-tau(z, zmax)/\mu)
        do j = ny, 1, -1
            p_tau(:, j) = exp(p_tau(:, j))
        end do
        !  p_tau = dexp(p_tau)         seg-fault; need ulimit -u unlimited

        ! emission function
        dummy = 2.0_wp*pi_wp*mu
        bcs = p_ab(:, 1)*dummy          ! save for calculation of surface flux
        p_ab = p_a*p_ab*dummy           ! absorption coefficient times emission function

! ###################################################################
        ! downward flux
        fd = -infrared%parameters(1)
        ! do j = ny, 1, -1
        !     tmp2(:, j) = p_ab(:, j)/p_tau(:, j)*p_tau(:, 1)
        ! end do
        tmp2 = p_ab/p_tau
        p_flux(:, ny) = 0.0_wp                                  ! boundary condition
        call OPR_Integral1(nxz, g, tmp2, p_flux, BCS_MAX)       ! recall this gives the negative of the integral
!        p_flux = fd*p_tau + p_flux                              ! negative going down
        ! do j = ny, 1, -1
        !     p_flux(:, j) = p_tau(:, j)*(fd + p_flux(:, j)/p_tau(:, 1))
        ! end do
        p_flux = p_tau*(fd + p_flux)

        ! calculate upward flux at the surface
        epsilon = infrared%parameters(4)
        bcs = epsilon*bcs - (1.0_wp - epsilon)*p_flux(:, 1)

        ! upward flux
        if (present(flux)) then
            ! tmp2 = p_ab*p_tau
            do j = ny, 1, -1
                tmp2(:, j) = p_ab(:, j)*p_tau(:, j)/p_tau(:, 1)
            end do
            flux(:, ny) = 0.0_wp                                ! boundary condition
            call OPR_Integral1(nxz, g, tmp2, flux, BCS_MAX)     ! recall this gives the negative of the integral
            do j = ny, 1, -1
                p_source(:, j) = p_a(:, j)*(p_tau(:, 1)/p_tau(:, j)*(bcs(:) + flux(:, j) - flux(:, 1)) - p_flux(:, j))- 2.0_wp*p_ab(:,j)
                ! p_flux(:, j) = bcs(:)*p_tau(:, 1)/p_tau(:, j) + flux(:, j) - flux(:, 1) &
                !                + p_flux(:, j)
                p_flux(:, j) = p_tau(:, 1)/p_tau(:, j)*(bcs(:) + flux(:, j) - flux(:, 1)) + p_flux(:, j)
                ! p_flux(:, j) = bcs(:)*p_tau(:, 1)/p_tau(:, j) + flux(:, j) - flux(:, 1) ! only up, testing
                ! p_flux(:, j) = p_tau(:, 1)/p_tau(:, j)*(bcs(:) + flux(:, j) - flux(:, 1))! only up, testing
                ! p_flux(:, j) = -p_flux(:, j) ! only down, testing
                p_source(:, j) = p_a(:, j)
            end do

#ifdef USE_ESSL
            call DGETMO(p_flux, nz, nz, nxy, flux, nxy)
#else
            call DNS_TRANSPOSE(p_flux, nz, nxy, nz, flux, nxy)
#endif
        end if

! ###################################################################
! Calculate heating rate
        fd = infrared%parameters(1)
        ! do j = ny, 1, -1
        !     p_source(:, j) = (fd*p_a(:, j) - p_ab(:, j))*p_tau(:, j) &
        !                      + (bcs(:)*p_a(:, j) - p_ab(:, j))*p_tau(:, 1)/p_tau(:, j)
        ! end do

#ifdef USE_ESSL
        call DGETMO(p_source, nz, nz, nxy, a_source, nxy)
#else
        call DNS_TRANSPOSE(p_source, nz, nxy, nz, a_source, nxy)
#endif

! ###################################################################
        nullify (p_a, p_tau, p_ab, p_source, p_flux)

        return
    end subroutine IR_Bulk1D

!########################################################################
!########################################################################
    subroutine Radiation_Infrared(infrared, nx, ny, nz, g, s, source, b, tmp1, tmp2, flux)
        use THERMO_ANELASTIC
        type(term_dt), intent(in) :: infrared
        integer(wi), intent(in) :: nx, ny, nz
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: s(nx*ny*nz, inb_scal_array)
        real(wp), intent(out) :: source(nx*ny*nz)   ! also used for bulk absorption coefficient
        real(wp), intent(inout) :: b(nx*ny*nz)      ! emission function
        real(wp), intent(inout) :: tmp1(nx*ny*nz), tmp2(nx*ny*nz)
        real(wp), intent(out), optional :: flux(nx*ny*nz)

        ! -----------------------------------------------------------------------
        real(wp) kappa, kappal, kappav

        !########################################################################
        select case (infrared%type)
        case (TYPE_IR_BULK1D_LIQUID, TYPE_BULK1DLOCAL)
            ! bulk absorption coefficient; we save it in array source to save memory
            kappa = infrared%parameters(2)
            source = kappa*s(:, infrared%scalar(1))
            if (imode_eqns == DNS_EQNS_ANELASTIC) then
                call THERMO_ANELASTIC_WEIGHT_INPLACE(nx, ny, nz, rbackground, source)
            end if

            if (present(flux)) then
                call IR_Bulk1D_Liquid(infrared, nx, ny, nz, g, source, flux)
            else
                call IR_Bulk1D_Liquid(infrared, nx, ny, nz, g, source)
            end if

        case (TYPE_IR_BULK1D)
            ! bulk absorption coefficient; we save it in array source to save memory
            kappal = infrared%parameters(2)
            kappav = infrared%parameters(3)
            source = kappal*s(:, infrared%scalar(1)) + kappav*(s(:, 2) - s(:, infrared%scalar(1)))
            if (imode_eqns == DNS_EQNS_ANELASTIC) then
                call THERMO_ANELASTIC_WEIGHT_INPLACE(nx, ny, nz, rbackground, source)
                ! calculate temperature to calculate emission function
                call THERMO_ANELASTIC_TEMPERATURE(nx, ny, nz, s, wrk3d)
            else
                ! tbd
            end if
            ! emission function
            b = sigma_o_pi*wrk3d**4.0_wp

            if (present(flux)) then
                call IR_Bulk1D(infrared, nx, ny, nz, g, source, b, tmp1, tmp2, flux)
            else
                call IR_Bulk1D(infrared, nx, ny, nz, g, source, b, tmp1, tmp2)
            end if

        end select

    end subroutine Radiation_Infrared

end module
