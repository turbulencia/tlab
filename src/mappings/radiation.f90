#include "dns_const.h"
#include "dns_error.h"

module Radiation
    use TLAB_CONSTANTS, only: wp, wi, pi_wp, BCS_MAX, BCS_MIN, efile, MAX_PROF
    use TLAB_TYPES, only: term_dt, grid_dt
    use TLAB_VARS, only: imode_eqns, inb_scal_array, isize_field
    use TLAB_ARRAYS, only: wrk2d, wrk3d
    use TLAB_POINTERS_3D, only: p_wrk3d
    use TLAB_PROCS, only: TLAB_WRITE_ASCII, TLAB_STOP, TLAB_ALLOCATE_ARRAY_DOUBLE
    use Thermodynamics, only: imixture
    use OPR_ODES
    use Integration
    implicit none
    private

    integer, parameter :: TYPE_NONE = 0
    integer, parameter :: TYPE_IR_GRAY_LIQUID = 1
    integer, parameter :: TYPE_IR_GRAY = 2
    integer, parameter :: TYPE_IR_3BANDS = 3
    integer, parameter :: TYPE_BULK1DLOCAL = 10         ! backwards compatibility, to be removed

    real(wp), parameter :: sigma = 5.67037442e-8_wp     ! Stefan-Boltzmann constant, W /m^2 /K^4
    real(wp) :: mu                                      ! mean direction parameter
    real(wp) :: epsilon                                 ! surface emissivity at ymin
    integer, parameter :: nbands_max = 3                ! maximum number of spectral bands
    integer :: nbands                                   ! number of spectral bands
    real(wp) beta(3, nbands_max)                        ! polynomial coefficients for band functions; assuming second-order polynomial
    real(wp) kappal(nbands_max), kappav(nbands_max)     ! mass absorption coefficients for liquid and vapor, for clarity
    real(wp), allocatable, target :: bcs_ht(:)          ! flux boundary condition at the top of the domain
    real(wp), allocatable, target :: bcs_hb(:)          ! flux boundary condition at the bottom of the domain
    real(wp), allocatable, target :: t_ht(:)            ! temperature at the top of the domain
    real(wp), allocatable, target :: tmp_rad(:, :)      ! 3D temporary arrays for radiation routine

    real(wp), pointer :: p_tau(:, :) => null()

    public :: Radiation_Initialize
    public :: Radiation_Infrared_Y

contains
!########################################################################
!########################################################################
    subroutine Radiation_Initialize(inifile)
        use TLAB_VARS, only: infrared, imax, kmax
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=512) sRes
        integer(wi) idummy, iband
        integer(wi) :: inb_tmp_rad = 0

        !########################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Infrared'

        call TLAB_WRITE_ASCII(bakfile, '#')
        call TLAB_WRITE_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLAB_WRITE_ASCII(bakfile, '#Type=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Scalar=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

        call SCANINICHAR(bakfile, inifile, block, 'Type', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') &
            call SCANINICHAR(bakfile, inifile, 'Main', 'TermRadiation', 'None', sRes)               ! backwards compatibility, to be removed
        if (trim(adjustl(sRes)) == 'none') then; infrared%type = TYPE_NONE
        else if (trim(adjustl(sRes)) == 'grayliquid') then; infrared%type = TYPE_IR_GRAY_LIQUID
        else if (trim(adjustl(sRes)) == 'gray') then; infrared%type = TYPE_IR_GRAY
        else if (trim(adjustl(sRes)) == 'threebands') then; infrared%type = TYPE_IR_3BANDS
        else if (trim(adjustl(sRes)) == 'bulk1dlocal') then; infrared%type = TYPE_BULK1DLOCAL    ! backwards compatibility, to be removed
        else
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Error in Radiation.Type.')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if

        infrared%active = .false.
        if (infrared%type /= TYPE_NONE) then
            call SCANINIINT(bakfile, inifile, block, 'Scalar', '1', idummy)
            infrared%active(idummy) = .true.

            infrared%parameters(:) = 0.0_wp
            call SCANINICHAR(bakfile, inifile, block, 'Parameters', '1.0', sRes)
            idummy = MAX_PROF
            call LIST_REAL(sRes, idummy, infrared%parameters)

        end if

        ! -------------------------------------------------------------------
        infrared%scalar = inb_scal_array                        ! By default, radiation is caused by last scalar

        if (imixture == MIXT_TYPE_AIRWATER .or. imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            if (infrared%type /= EQNS_NONE) then
                infrared%active(inb_scal_array) = .true.        ! liquid
                infrared%active(inb_scal_array + 1) = .true.    ! buoyancy
            end if

        end if

        if (infrared%type == TYPE_BULK1DLOCAL) then             ! backwards compatibility
            infrared%parameters(1) = infrared%parameters(1)*infrared%parameters(2)
            infrared%parameters(3) = infrared%parameters(3)*infrared%parameters(2)
            infrared%parameters(2) = 1.0_wp/infrared%parameters(2)
            infrared%type = TYPE_IR_GRAY_LIQUID
        end if

        select case (infrared%type)
        case (TYPE_IR_GRAY_LIQUID)
            nbands = 1

            kappal(1) = infrared%parameters(2)     ! mass absorption coefficient of liquid
            ! infrared%parameters(3) upward flux at domain bottom

        case (TYPE_IR_GRAY)
            nbands = 1

            kappal(1) = infrared%parameters(2)      ! mass absorption coefficient of liquid
            kappav(1) = infrared%parameters(3)      ! mass absorption coefficient of vapor
            epsilon = infrared%parameters(4)        ! surface emissivity at ymin

        case (TYPE_IR_3BANDS)
            nbands = 3

            ! For the airwater mixture
            kappal(1:3) = infrared%parameters(2)    ! mass absorption coefficient of liquid, same in all bands
            kappav(1) = infrared%parameters(3)      ! mass absorption coefficient of vapor, band 1
            kappav(2) = infrared%parameters(4)      ! mass absorption coefficient of vapor, band 2
            kappav(3) = 0.0_wp                      ! assume band 3 is defined by vapor being transparent
            epsilon = infrared%parameters(5)        ! surface emissivity at ymin

            beta(1:3, 1) = [2.6774e-1_wp, -1.3344e-3_wp, 1.8017e-6_wp] ! coefficients for band 1
            beta(1:3, 2) = [-2.2993e-2_wp, 8.7439e-5_wp, 1.4744e-7_wp] ! coefficients for band 2
            beta(1:3, nbands) = [1.0_wp, 0.0_wp, 0.0_wp]               ! last band from equation sum beta_i = 1
            do iband = 1, nbands - 1
                beta(1:3, nbands) = beta(1:3, nbands) - beta(1:3, iband)
            end do

            inb_tmp_rad = 5                         ! Additional memory space

        end select

        mu = 0.5_wp*(1.0_wp/sqrt(3.0_wp) + 1.0_wp/sqrt(2.0_wp))     ! mean direction, in (1/sqrt{3}, 1/sqrt{2})
        ! mu = 1.0_wp/sqrt(2.0_wp)
        ! mu = 0.5_wp     ! testing

        ! !########################################################################
        ! bakfile = trim(adjustl(inifile))//'.bak'
        ! block = 'Visible'

        !########################################################################
        ! Local allocation
        allocate (bcs_ht(imax*kmax), bcs_hb(imax*kmax), t_ht(imax*kmax))
        call TLAB_ALLOCATE_ARRAY_DOUBLE(__FILE__, tmp_rad, [isize_field, inb_tmp_rad], 'tmp-rad')

        return
    end subroutine Radiation_Initialize

!########################################################################
!########################################################################
    subroutine Radiation_Infrared_Y(infrared, nx, ny, nz, g, s, source, b, tmp1, tmp2, flux)
        use THERMO_ANELASTIC
        type(term_dt), intent(in) :: infrared
        integer(wi), intent(in) :: nx, ny, nz
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: s(nx*ny*nz, inb_scal_array)
        real(wp), intent(out) :: source(nx*ny*nz)           ! also used for absorption coefficient
        real(wp), intent(inout) :: b(nx*ny*nz)              ! emission function
        real(wp), intent(inout) :: tmp1(nx*ny*nz), tmp2(nx*ny*nz)
        real(wp), intent(out), optional :: flux(nx*ny*nz)

        target :: source, b, tmp1, flux

        ! -----------------------------------------------------------------------
        integer iband, nxy, nxz
        real(wp), dimension(:, :), pointer :: p_bcs
        real(wp), pointer :: p_source(:) => null()
        real(wp), pointer :: p_b(:) => null()
        real(wp), pointer :: p_flux_down(:) => null()
        real(wp), pointer :: p_flux_up(:) => null()

        !########################################################################
        nxy = nx*ny     ! For transposition to make y direction the last one
        nxz = nx*nz

        p_source => b
        p_b => source
        p_flux_down => tmp1
        if (present(flux)) then
            p_flux_up => flux
        end if

        p_tau(1:nxz, 1:ny) => wrk3d(1:nxz*ny)                   ! set pointer to optical depth and transmission functions used in routines below

        ! -----------------------------------------------------------------------
        select case (infrared%type)
        case (TYPE_IR_GRAY_LIQUID)
            wrk3d = kappal(1)*s(:, infrared%scalar(1))          ! absorption coefficient in array source to save memory
            if (imode_eqns == DNS_EQNS_ANELASTIC) then
                call THERMO_ANELASTIC_WEIGHT_INPLACE(nx, ny, nz, rbackground, wrk3d)
            end if
            ! Local transposition: make x-direction the last one. Same in the similar blocks below
#ifdef USE_ESSL
            call DGETMO(wrk3d, nxy, nxy, nz, p_source, nz)
#else
            call DNS_TRANSPOSE(wrk3d, nxy, nz, nxy, p_source, nz)
#endif

            bcs_ht = infrared%parameters(1)                     ! downward flux at domain top
            if (present(flux)) then                             ! solve radiative transfer equation along y
                call IR_RTE1_Liquid(infrared, nxz, ny, g, p_source, p_flux_down, p_flux_up)
            else
                call IR_RTE1_Liquid(infrared, nxz, ny, g, p_source)
            end if

            ! -----------------------------------------------------------------------
        case (TYPE_IR_GRAY)
            if (imode_eqns == DNS_EQNS_ANELASTIC) then
                call THERMO_ANELASTIC_TEMPERATURE(nx, ny, nz, s, wrk3d)
            else
                ! tbd
            end if
            wrk3d = sigma*wrk3d**4.0_wp                         ! emission function, Stefan-Boltzmann law
#ifdef USE_ESSL
            call DGETMO(wrk3d, nxy, nxy, nz, p_b, nz)
#else
            call DNS_TRANSPOSE(wrk3d, nxy, nz, nxy, p_b, nz)
#endif

            wrk3d = kappal(1)*s(:, infrared%scalar(1)) + kappav(1)*(s(:, 2) - s(:, infrared%scalar(1))) ! absorption coefficient
            if (imode_eqns == DNS_EQNS_ANELASTIC) then
                call THERMO_ANELASTIC_WEIGHT_INPLACE(nx, ny, nz, rbackground, wrk3d)
            end if
#ifdef USE_ESSL
            call DGETMO(wrk3d, nxy, nxy, nz, p_source, nz)
#else
            call DNS_TRANSPOSE(wrk3d, nxy, nz, nxy, p_source, nz)
#endif

            bcs_ht = infrared%parameters(1)                     ! downward flux at domain top

            if (present(flux)) then                             ! solve radiative transfer equation along y
                call IR_RTE1_Global(infrared, nxz, ny, g, p_source, p_b, p_flux_down, p_flux_up)
                ! call IR_RTE1_Local(infrared, nxz, ny, g, p_source, p_b, p_flux_down, tmp2, p_flux_up)
                ! call IR_RTE1_Incremental(infrared, nxz, ny, g, p_source, p_b, p_flux_down, p_flux_up)
            else
                call IR_RTE1_Global(infrared, nxz, ny, g, p_source, p_b, p_flux_down, tmp2)
                ! call IR_RTE1_Local(infrared, nxz, ny, g, p_source, p_b, p_flux_down, tmp2)
                ! call IR_RTE1_Incremental(infrared, nxz, ny, g, p_source, p_b, p_flux_down)
            end if

            ! -----------------------------------------------------------------------
        case (TYPE_IR_3BANDS)
            if (imode_eqns == DNS_EQNS_ANELASTIC) then
                call THERMO_ANELASTIC_TEMPERATURE(nx, ny, nz, s, wrk3d)
            else
                ! tbd
            end if
            p_wrk3d(1:nx, 1:ny, 1:nz) => tmp_rad(1:nx*ny*nz, 1) ! save T at top boundary
            p_bcs(1:nx, 1:nz) => t_ht(1:nx*nz)
            p_bcs(1:nx, 1:nz) = p_wrk3d(1:nx, ny, 1:nz)
#ifdef USE_ESSL
            call DGETMO(wrk3d, nxy, nxy, nz, tmp_rad, nz)
#else
            call DNS_TRANSPOSE(wrk3d, nxy, nz, nxy, tmp_rad, nz)
#endif

            tmp_rad(:, 2) = s(:, 2) - s(:, infrared%scalar(1))  ! vapor

            ! last band
            iband = nbands
            p_b = sigma*tmp_rad(:, 1)**4.0_wp*(beta(1, iband) + tmp_rad(:, 1)*(beta(2, iband) + tmp_rad(:, 1)*beta(3, iband)))

            wrk3d = kappal(iband)*s(:, infrared%scalar(1))
            if (imode_eqns == DNS_EQNS_ANELASTIC) then
                call THERMO_ANELASTIC_WEIGHT_INPLACE(nx, ny, nz, rbackground, wrk3d)
            end if
#ifdef USE_ESSL
            call DGETMO(wrk3d, nxy, nxy, nz, p_source, nz)
#else
            call DNS_TRANSPOSE(wrk3d, nxy, nz, nxy, p_source, nz)
#endif

            bcs_ht = infrared%parameters(1)*(beta(1, iband) + t_ht*(beta(2, iband) + t_ht*beta(3, iband)))

            if (present(flux)) then
                call IR_RTE1_Global(infrared, nxz, ny, g, p_source, p_b, p_flux_down, p_flux_up)
                ! call IR_RTE1_Local(infrared, nxz, ny, g, p_source, p_b, p_flux_down, tmp2, p_flux_up)
                ! call IR_RTE1_Incremental(infrared, nxz, ny, g, p_source, p_b, p_flux_down, p_flux_up)
            else
                call IR_RTE1_Global(infrared, nxz, ny, g, p_source, p_b, p_flux_down, tmp2)
                ! call IR_RTE1_Local(infrared, nxz, ny, g, p_source, p_b, p_flux_down, tmp2)
                ! call IR_RTE1_Incremental(infrared, nxz, ny, g, p_source, p_b, p_flux_down)
            end if

            ! the first nband-1 bands
            do iband = 1, nbands - 1
                tmp_rad(:, 5) = sigma*tmp_rad(:, 1)**4.0_wp*(beta(1, iband) + tmp_rad(:, 1)*(beta(2, iband) + tmp_rad(:, 1)*beta(3, iband)))

                wrk3d = kappal(iband)*s(:, infrared%scalar(1)) + kappav(iband)*tmp_rad(:, 2)
                if (imode_eqns == DNS_EQNS_ANELASTIC) then
                    call THERMO_ANELASTIC_WEIGHT_INPLACE(nx, ny, nz, rbackground, wrk3d)        ! multiply by density
                end if
#ifdef USE_ESSL
                call DGETMO(wrk3d, nxy, nxy, nz, tmp_rad(:, 3), nz)
#else
                call DNS_TRANSPOSE(wrk3d, nxy, nz, nxy, tmp_rad(:, 3), nz)
#endif

                bcs_ht = infrared%parameters(1)*(beta(1, iband) + t_ht*(beta(2, iband) + t_ht*beta(3, iband)))  ! downward flux at domain top

                if (present(flux)) then             ! solve radiative transfer equation along y
                    call IR_RTE1_Global(infrared, nxz, ny, g, tmp_rad(:, 3), tmp_rad(:, 5), tmp_rad(:, 4), p_b)
                    ! call IR_RTE1_Local(infrared, nxz, ny, g, tmp_rad(:, 3), tmp_rad(:, 5), tmp_rad(:, 4), tmp2, p_b)
                    ! call IR_RTE1_Incremental(infrared, nxz, ny, g, tmp_rad(:, 3), tmp_rad(:, 5), tmp_rad(:, 4), p_b)
                    p_flux_down = p_flux_down + tmp_rad(:, 4)
                    p_flux_up = p_flux_up + p_b
                else
                    call IR_RTE1_Global(infrared, nxz, ny, g, tmp_rad(:, 3), tmp_rad(:, 5), tmp_rad(:, 4), tmp2)
                    ! call IR_RTE1_Local(infrared, nxz, ny, g, tmp_rad(:, 3), tmp_rad(:, 5), tmp_rad(:, 4), tmp2)
                    ! call IR_RTE1_Incremental(infrared, nxz, ny, g, tmp_rad(:, 3), tmp_rad(:, 5), tmp_rad(:, 4))
                end if

                p_source = p_source + tmp_rad(:, 3)

            end do

        end select

        !########################################################################
        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(p_source, nz, nz, nxy, source, nxy)
        if (present(flux)) then
            p_flux_down = p_flux_up - p_flux_down                 ! net flux upwards = upward flux - downward flux
            call DGETMO(p_flux_up, nz, nz, nxy, b, nxy)
            call DGETMO(p_flux_down, nz, nz, nxy, flux, nxy)
        end if
#else
        call DNS_TRANSPOSE(p_source, nz, nxy, nz, source, nxy)
        if (present(flux)) then
            p_flux_down = p_flux_up - p_flux_down
            call DNS_TRANSPOSE(p_flux_up, nz, nxy, nz, b, nxy)
            call DNS_TRANSPOSE(p_flux_down, nz, nxy, nz, flux, nxy)
        end if
#endif

        nullify (p_source, p_b, p_flux_down)
        if (associated(p_flux_up)) nullify (p_flux_up)
        nullify (p_tau)

    end subroutine Radiation_Infrared_Y

!########################################################################
! Solve radiative transfer equation along 1 direction
! We do not treat separately 2d and 3d cases for the transposition because it was a bit messy...
!########################################################################
    ! only liquid
    subroutine IR_RTE1_Liquid(infrared, nlines, ny, g, a_source, flux_down, flux_up)
        type(term_dt), intent(in) :: infrared
        integer(wi), intent(in) :: nlines, ny
        type(grid_dt), intent(in) :: g
        real(wp), intent(inout) :: a_source(nlines, ny)      ! input as bulk absorption coefficent, output as source
        real(wp), intent(out), optional :: flux_down(nlines, ny), flux_up(nlines, ny)

! -----------------------------------------------------------------------
        integer(wi) j

! #######################################################################
        ! calculate f_j = exp(-tau(z, zmax)/\mu)
        p_tau(:, ny) = 0.0_wp                                   ! boundary condition
        call OPR_Integral1(nlines, g, a_source, p_tau, BCS_MAX)         ! recall this gives the negative of the integral
        ! call Int_Trapezoidal_f(a_source, g%nodes, p_tau, BCS_MAX)
        ! call Int_Simpson_Biased_f(a_source, g%nodes, p_tau, BCS_MAX)
        do j = ny, 1, -1
            p_tau(:, j) = exp(p_tau(:, j))
        end do
        !  p_tau = dexp(p_tau)         seg-fault; need ulimit -u unlimited

        ! Calculate heating rate
        bcs_hb = infrared%parameters(3)
        if (abs(infrared%parameters(3)) > 0.0_wp) then
            do j = ny, 1, -1
                a_source(:, j) = a_source(:, j)*(p_tau(:, j)*bcs_ht(1:nlines) &                   ! downward flux
                                                 + p_tau(:, 1)/p_tau(:, j)*bcs_hb(1:nlines))      ! upward flux
            end do
        else
            do j = ny, 1, -1
                a_source(:, j) = a_source(:, j)*p_tau(:, j)*bcs_ht(1:nlines)
            end do
        end if

        ! Calculate flux, if necessary
        if (present(flux_up)) then
            do j = ny, 1, -1
                flux_down(:, j) = bcs_ht(1:nlines)*p_tau(:, j)                       ! downward flux
                flux_up(:, j) = bcs_hb(1:nlines)*p_tau(:, 1)/p_tau(:, j)        ! upward flux
            end do

        end if

        return
    end subroutine IR_RTE1_Liquid

!########################################################################
!########################################################################
    subroutine IR_RTE1_Incremental(infrared, nlines, ny, g, a_source, b, flux_down, flux_up)
        type(term_dt), intent(in) :: infrared
        integer(wi), intent(in) :: nlines, ny
        type(grid_dt), intent(in) :: g
        real(wp), intent(inout) :: a_source(nlines, ny)         ! input as bulk absorption coefficent, output as source
        real(wp), intent(inout) :: b(nlines, ny)                ! input as emission function, output as upward flux, if flux is to be return
        real(wp), intent(inout) :: flux_down(nlines, ny)             ! flux_down for intermediate calculations and net flux as output
        real(wp), intent(out), optional :: flux_up(nlines, ny)

! -----------------------------------------------------------------------
        integer(wi) j
        real(wp) dummy
        real(wp), pointer :: p_wrk2d_1(:) => null()
        real(wp), pointer :: p_wrk2d_2(:) => null()

! #######################################################################
        p_wrk2d_1(1:nlines) => wrk2d(1:nlines, 1)
        p_wrk2d_2(1:nlines) => wrk2d(1:nlines, 2)

! ###################################################################
        ! absorption coefficient; divide by mean direction
        dummy = 1.0_wp/mu
        a_source = a_source*dummy

        ! emission function
        bcs_hb(1:nlines) = b(1:nlines, 1)   ! save for calculation of surface flux
        b = b*a_source                      ! absorption coefficient times emission function

        ! transmission function I_{j-1,j}  = exp(-tau(z_{j-1}, z_j)/\mu)
        p_tau(:, 1) = 0.0_wp                                    ! boundary condition
        ! call OPR_Integral1(nxz, g, a_source, p_tau, BCS_MIN)
        ! call Int_Trapezoidal_f(a_source, g%nodes, p_tau, BCS_MIN)
        call Int_Simpson_Biased_f(a_source, g%nodes, p_tau, BCS_MIN)
        do j = ny, 2, -1
            p_tau(:, j) = exp(p_tau(:, j - 1) - p_tau(:, j))
        end do
        p_tau(:, 1) = 1.0_wp        ! this value is not used

        ! ###################################################################
        ! downward flux; positive going down
        j = ny
        flux_down(1:nlines, j) = bcs_ht(1:nlines)

        do j = ny - 1, 1, -1
            ! Integral contribution from emission function using a trapezoidal rule
            p_wrk2d_2 = b(:, j + 1)
            p_wrk2d_1 = b(:, j)/p_tau(:, j + 1)
            p_wrk2d_1 = 0.5_wp*(p_wrk2d_1 + p_wrk2d_2)*(g%nodes(j + 1) - g%nodes(j))

            flux_down(:, j) = p_tau(:, j + 1)*(flux_down(:, j + 1) + p_wrk2d_1)
        end do

        ! ###################################################################
        ! bottom boundary condition; calculate upward flux at the bottom
        bcs_hb(1:nlines) = epsilon*bcs_hb(1:nlines) + (1.0_wp - epsilon)*flux_down(:, 1)

        ! ###################################################################
        ! upward flux and source
        j = 1
        a_source(:, j) = a_source(:, j)*(bcs_hb(1:nlines) + flux_down(:, j)) - 2.0_wp*b(:, j)

        if (present(flux_up)) then                                          ! Save fluxes
            flux_up(:, j) = bcs_hb(1:nlines)                                ! upward flux

            do j = 2, ny
                ! Integral contribution from emission function using a trapezoidal rule
                p_wrk2d_1 = b(:, j - 1)
                p_wrk2d_2 = b(:, j)/p_tau(:, j)
                p_wrk2d_1 = 0.5_wp*(p_wrk2d_1 + p_wrk2d_2)*(g%nodes(j) - g%nodes(j - 1))

                bcs_hb = p_tau(:, j)*(bcs_hb(1:nlines) + p_wrk2d_1)
                a_source(:, j) = a_source(:, j)*(bcs_hb(1:nlines) + flux_down(:, j)) - 2.0_wp*b(:, j)

                flux_up(:, j) = bcs_hb(1:nlines)
            end do

        else                                                                ! Do not save fluxes
            do j = 2, ny
                ! Integral contribution from emission function using a trapezoidal rule
                p_wrk2d_1 = b(:, j - 1)
                p_wrk2d_2 = b(:, j)/p_tau(:, j)
                p_wrk2d_1 = 0.5_wp*(p_wrk2d_1 + p_wrk2d_2)*(g%nodes(j) - g%nodes(j - 1))

                bcs_hb = p_tau(:, j)*(bcs_hb + p_wrk2d_1)
                a_source(:, j) = a_source(:, j)*(bcs_hb + flux_down(:, j)) - 2.0_wp*b(:, j)
            end do

        end if

! ###################################################################
        nullify (p_wrk2d_1, p_wrk2d_2)

        return
    end subroutine IR_RTE1_Incremental

!########################################################################
!########################################################################
    subroutine IR_RTE1_Local(infrared, nlines, ny, g, a_source, b, flux_down, tmp2, flux_up)
        type(term_dt), intent(in) :: infrared
        integer(wi), intent(in) :: nlines, ny
        type(grid_dt), intent(in) :: g
        real(wp), intent(inout) :: a_source(nlines, ny)         ! input as bulk absorption coefficent, output as source
        real(wp), intent(inout) :: b(nlines, ny)                ! input as emission function, output as upward flux, if flux is to be return
        real(wp), intent(inout) :: flux_down(nlines, ny)             ! flux_down for intermediate calculations and net flux as output
        real(wp), intent(inout) :: tmp2(nlines, ny)
        real(wp), intent(out), optional :: flux_up(nlines, ny)

! -----------------------------------------------------------------------
        integer(wi) j, k
        real(wp) dummy
        real(wp), pointer :: p_flux_1(:) => null()
        real(wp), pointer :: p_flux_2(:) => null()

! #######################################################################
        p_flux_1(1:nlines) => wrk2d(1:nlines, 1)
        p_flux_2(1:nlines) => wrk2d(1:nlines, 2)

! ###################################################################
        ! absorption coefficient; divide by mean direction
        dummy = 1.0_wp/mu
        a_source = a_source*dummy

        ! emission function
        bcs_hb(1:nlines) = b(1:nlines, 1)   ! save for calculation of surface flux
        b = b*a_source                      ! absorption coefficient times emission function

        ! transmission function I_{j-1,j} = exp(-tau(z_{j-1}, z_j)/\mu)
        p_tau(:, 1) = 0.0_wp                                    ! boundary condition
        ! call OPR_Integral1(nxz, g, a_source, p_tau, BCS_MIN)
        ! call Int_Trapezoidal_f(a_source, g%nodes, p_tau, BCS_MIN)
        call Int_Simpson_Biased_f(a_source, g%nodes, p_tau, BCS_MIN)
        do j = ny, 2, -1
            p_tau(:, j) = exp(p_tau(:, j - 1) - p_tau(:, j))
        end do
        p_tau(:, 1) = 1.0_wp        ! this value is not used

        ! ###################################################################
        ! downward flux; positive going down
        j = ny
        p_flux_1(1:nlines) = bcs_ht(1:nlines)
        flux_down(:, j) = p_flux_1

        do j = ny - 1, 1, -1
            p_flux_1 = p_flux_1*p_tau(:, j + 1)                     ! accumulate transmission

            p_flux_2 = 1.0_wp                                       ! calculate emission; p_flux_2 is aux array in the next loop
            tmp2(:, j) = b(:, j)
            do k = j + 1, ny
                p_flux_2 = p_flux_2*p_tau(:, k)
                tmp2(:, k) = b(:, k)*p_flux_2
            end do
            call Int_Simpson_v(tmp2(:, j:ny), g%nodes(j:ny), flux_down(:, j))
            ! call Int_Trapezoidal_v(tmp2(:, j:ny), g%nodes(j:ny), flux_down(:, j))
            flux_down(:, j) = flux_down(:, j) + p_flux_1
        end do

        ! ###################################################################
        ! bottom boundary condition; calculate upward flux at the bottom
        bcs_hb(1:nlines) = epsilon*bcs_hb(1:nlines) + (1.0_wp - epsilon)*flux_down(:, 1)

        ! ###################################################################
        ! upward flux and net terms
        j = 1
        p_flux_1(1:nlines) = bcs_hb(1:nlines)
        a_source(:, j) = a_source(:, j)*(p_flux_1 + flux_down(:, j)) - 2.0_wp*b(:, j)

        if (present(flux_up)) then                                  ! I need additional space to store fluxes
            flux_up(:, j) = p_flux_1                                ! upward flux

            do j = 2, ny
                p_flux_1 = p_flux_1*p_tau(:, j)

                p_flux_2 = 1.0_wp       ! used as auxiliary array in the next loop
                tmp2(:, j) = b(:, j)
                do k = j - 1, 1, -1
                    p_flux_2 = p_flux_2*p_tau(:, k + 1)
                    tmp2(:, k) = b(:, k)*p_flux_2
                end do
                call Int_Simpson_v(tmp2(:, 1:j), g%nodes(1:j), p_flux_2)
                ! call Int_Trapezoidal_v(tmp2(:, 1:j), g%nodes(1:j), p_flux_2)
                a_source(:, j) = a_source(:, j)*(p_flux_2 + p_flux_1 + flux_down(:, j)) - 2.0_wp*b(:, j)

                flux_up(:, j) = p_flux_2 + p_flux_1                 ! upward flux
            end do

        else                ! we only need the source term
            do j = 2, ny
                p_flux_1 = p_flux_1*p_tau(:, j)                     ! accumulate transmission

                p_flux_2 = 1.0_wp                                   ! calculate emission; p_flux_2 is aux array in the next loop
                tmp2(:, j) = b(:, j)
                do k = j - 1, 1, -1
                    p_flux_2 = p_flux_2*p_tau(:, k + 1)
                    tmp2(:, k) = b(:, k)*p_flux_2
                end do
                call Int_Simpson_v(tmp2(:, 1:j), g%nodes(1:j), p_flux_2)
                ! call Int_Trapezoidal_v(tmp2(:, 1:j), g%nodes(1:j), p_flux_2)
                a_source(:, j) = a_source(:, j)*(p_flux_2 + p_flux_1 + flux_down(:, j)) - 2.0_wp*b(:, j)

            end do

        end if

! #######################################################################
        nullify (p_flux_1, p_flux_2)

        return
    end subroutine IR_RTE1_Local

!########################################################################
!########################################################################
    subroutine IR_RTE1_Global(infrared, nlines, ny, g, a_source, b, flux_down, flux_up)
        type(term_dt), intent(in) :: infrared
        integer(wi), intent(in) :: nlines, ny
        type(grid_dt), intent(in) :: g
        real(wp), intent(inout) :: a_source(nlines, ny)         ! input as bulk absorption coefficent, output as source
        real(wp), intent(inout) :: b(nlines, ny)                ! input as emission function, output as upward flux, if flux is to be return
        real(wp), intent(inout) :: flux_down(nlines, ny)             ! flux_down for intermediate calculations and net flux as output
        real(wp), intent(inout) :: flux_up(nlines, ny)

! -----------------------------------------------------------------------
        integer(wi) j
        real(wp) dummy

! ###################################################################
        ! absorption coefficient; divide by mean direction
        dummy = 1.0_wp/mu
        a_source = a_source*dummy

        ! emission function
        bcs_hb = b(:, 1)                ! save for calculation of surface flux
        b = b*a_source                  ! absorption coefficient times emission function

        ! ###################################################################
        ! downward flux; positive going down

        ! transmission function I_j = exp(-tau(z_j, zmax)/\mu)
        p_tau(:, ny) = 0.0_wp                                   ! boundary condition
        ! call OPR_Integral1(nlines, g, a_source, p_tau, BCS_MAX)         ! recall this gives the negative of the integral
        ! call Int_Trapezoidal_f(a_source, g%nodes, p_tau, BCS_MAX)
        call Int_Simpson_Biased_f(a_source, g%nodes, p_tau, BCS_MAX)
        do j = ny, 1, -1
            p_tau(:, j) = exp(-p_tau(:, j))
        end do
        !  p_tau = dexp(p_tau)         seg-fault; need ulimit -u unlimited

        flux_down = b/p_tau
        ! call Int_Trapezoidal_Increments_InPlace(flux_down, g%nodes, BCS_MAX)                   ! Calculate I_j = int_{x_{j}}^{x_{j+1}}
        call Int_Simpson_Biased_Increments_InPlace(flux_down, g%nodes, wrk2d(:, 1), BCS_MAX)   ! Calculate I_j = int_{x_{j}}^{x_{j+1}}
        j = ny
        wrk2d(1:nlines, 1) = 0.0_wp                                         ! accumulate emission; using wrk2d as aux array
        flux_down(:, j) = p_tau(:, j)*(bcs_ht(1:nlines) + wrk2d(1:nlines, 1))
        do j = ny - 1, 1, -1
            wrk2d(1:nlines, 1) = wrk2d(1:nlines, 1) + flux_down(:, j)            ! accumulate emission
            flux_down(:, j) = p_tau(:, j)*(bcs_ht(1:nlines) + wrk2d(1:nlines, 1))
        end do

        ! ###################################################################
        ! bottom boundary condition; calculate upward flux at the bottom
        bcs_hb = epsilon*bcs_hb + (1.0_wp - epsilon)*flux_down(:, 1)

        ! ###################################################################
        ! upward flux and net terms

        ! transmission function I_j = exp(-tau(zmin, z)/\mu)
        p_tau(:, 1) = 0.0_wp                                                ! boundary condition
        ! call OPR_Integral1(nlines, g, a_source, p_tau, BCS_MIN)
        ! call Int_Trapezoidal_f(a_source, g%nodes, p_tau, BCS_MIN)
        call Int_Simpson_Biased_f(a_source, g%nodes, p_tau, BCS_MIN)
        do j = 1, ny
            p_tau(:, j) = exp(-p_tau(:, j))
        end do

        flux_up = b/p_tau
        ! call Int_Trapezoidal_Increments_InPlace(flux_up, g%nodes, BCS_MIN)                     ! Calculate I_j = int_{x_{j-1}}^{x_{j}}
        call Int_Simpson_Biased_Increments_InPlace(flux_up, g%nodes, wrk2d(:, 1), BCS_MIN)     ! Calculate I_j = int_{x_{j-1}}^{x_{j}}
        j = 1
        wrk2d(1:nlines, 1) = 0.0_wp                                         ! accumulate emission; using wrk2d as aux array
        flux_up(:, j) = p_tau(:, j)*(bcs_hb(1:nlines) + wrk2d(1:nlines, 1))
        !
        a_source(:, j) = a_source(:, j)*(flux_up(:, j) + flux_down(:, j)) - 2.0_wp*b(:, j)
        do j = 2, ny
            wrk2d(1:nlines, 1) = wrk2d(1:nlines, 1) + flux_up(:, j)         ! accumulate emission
            flux_up(:, j) = p_tau(:, j)*(bcs_hb(1:nlines) + wrk2d(1:nlines, 1))
            !
            a_source(:, j) = a_source(:, j)*(flux_up(:, j) + flux_down(:, j)) - 2.0_wp*b(:, j)
        end do

        return
    end subroutine IR_RTE1_Global

end module
