#include "dns_const.h"
#include "dns_error.h"

! Compressible formulation uses simply the gravity force rho *g
! Incompressible formulation uses buoyancy and different forms of buoyancy function
! Anelastic formulation uses (rho-rho_ref) *g and scaleheight defines rho_ref

module Gravity
    use TLab_Constants, only: wp, wi, small_wp, efile, lfile, wfile, MAX_PROF
    use TLab_Types, only: term_dt
    use TLAB_VARS, only: inb_scal, inb_scal_array
    use TLAB_VARS, only: froude
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none
    private

    type(term_dt) :: buoyancy
    real(wp), allocatable :: bbackground(:)

    ! integer, parameter :: EQNS_BOD_HOMOGENEOUS = 5
    ! integer, parameter :: EQNS_BOD_LINEAR = 6
    ! integer, parameter :: EQNS_BOD_BILINEAR = 7
    ! integer, parameter :: EQNS_BOD_QUADRATIC = 8
    ! integer, parameter :: EQNS_BOD_NORMALIZEDMEAN = 9
    ! integer, parameter :: EQNS_BOD_SUBTRACTMEAN = 10

    public :: buoyancy, bbackground
    public :: Gravity_Initialize
    public :: Gravity_Hydrostatic_Enthalpy
    public :: Gravity_Buoyancy, Gravity_Buoyancy_Source

contains
    !########################################################################
    !########################################################################
    subroutine Gravity_Initialize(inifile)
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=512) sRes
        integer(wi) idummy

        !########################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Gravity'

        call ScanFile_Char(bakfile, inifile, block, 'Vector', 'void', sRes)                             ! backwards compatibility, to be removed
        if (trim(adjustl(sRes)) == 'void') block = 'BodyForce'

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<none/Explicit/Homogeneous/Linear/Bilinear/Quadratic>')
        call TLab_Write_ASCII(bakfile, '#Vector=<Gx,Gy,Gz>')
        call TLab_Write_ASCII(bakfile, '#Parameters=<value>')

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') &
            call ScanFile_Char(bakfile, inifile, 'Main', 'TermBodyForce', 'none', sRes)                 ! backwards compatibility, to be removed
        if (trim(adjustl(sRes)) == 'none') then; buoyancy%type = EQNS_NONE
        else if (trim(adjustl(sRes)) == 'explicit') then; buoyancy%type = EQNS_EXPLICIT
        else if (trim(adjustl(sRes)) == 'homogeneous') then; buoyancy%type = EQNS_BOD_HOMOGENEOUS
        else if (trim(adjustl(sRes)) == 'linear') then; buoyancy%type = EQNS_BOD_LINEAR
        else if (trim(adjustl(sRes)) == 'bilinear') then; buoyancy%type = EQNS_BOD_BILINEAR
        else if (trim(adjustl(sRes)) == 'quadratic') then; buoyancy%type = EQNS_BOD_QUADRATIC
        else if (trim(adjustl(sRes)) == 'normalizedmean') then; buoyancy%type = EQNS_BOD_NORMALIZEDMEAN
        else if (trim(adjustl(sRes)) == 'subtractmean') then; buoyancy%type = EQNS_BOD_SUBTRACTMEAN
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong TermBodyForce option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        if (any([EQNS_BOD_LINEAR, EQNS_BOD_BILINEAR, EQNS_BOD_QUADRATIC] == buoyancy%type) .and. inb_scal == 0) then
            call TLab_Write_ASCII(wfile, __FILE__//'. Zero scalars; setting TermBodyForce equal to none.')
            buoyancy%type = EQNS_NONE
        end if

        buoyancy%vector = 0.0_wp
        call ScanFile_Char(bakfile, inifile, block, 'Vector', '0.0,0.0,0.0', sRes)
        idummy = 3
        call LIST_REAL(sRes, idummy, buoyancy%vector)

        buoyancy%active = .false.
        if (abs(buoyancy%vector(1)) > 0.0_wp) then; buoyancy%active(1) = .true.; call TLab_Write_ASCII(lfile, 'Gravity along Ox.'); end if
        if (abs(buoyancy%vector(2)) > 0.0_wp) then; buoyancy%active(2) = .true.; call TLab_Write_ASCII(lfile, 'Gravity along Oy.'); end if
        if (abs(buoyancy%vector(3)) > 0.0_wp) then; buoyancy%active(3) = .true.; call TLab_Write_ASCII(lfile, 'Gravity along Oz.'); end if

        if (froude > 0.0_wp) then
            buoyancy%vector(:) = buoyancy%vector(:)/froude ! adding the froude number into the vector g
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Froude number must be nonzero if buoyancy is retained.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        if (buoyancy%type /= EQNS_NONE) then

            buoyancy%parameters(:) = 0.0_wp
            call ScanFile_Char(bakfile, inifile, block, 'Parameters', '0.0', sRes)
            idummy = MAX_PROF
            call LIST_REAL(sRes, idummy, buoyancy%parameters)
            buoyancy%scalar(1) = idummy

        end if

        return
    end subroutine Gravity_Initialize

    !########################################################################
    ! Compute hydrostatic equilibrium from profiles s=(h,q_t) where h is the enthalpy
    ! Evaluate the integral \int_yref^y dx/H(x), where H(x) is the scale height in the system
    !########################################################################
    subroutine Gravity_Hydrostatic_Enthalpy(g, s, ep, T, p, yref, pref, wrk1d)
        use TLab_Constants, only: BCS_MIN
        use FDM, only: grid_dt
        use TLAB_VARS, only: damkohler
        use Thermodynamics
        use THERMO_ANELASTIC
        use THERMO_AIRWATER
        use THERMO_THERMAL
        use OPR_ODES

        type(grid_dt), intent(in) :: g
        real(wp), dimension(g%size, inb_scal_array), intent(inout) :: s      ! We calculate equilibrium composition
        real(wp), dimension(g%size), intent(out) :: ep, T, p
        real(wp), intent(in) :: yref, pref
        real(wp), dimension(g%size, 2), intent(inout) :: wrk1d

        ! -------------------------------------------------------------------
        integer(wi) iter, niter, j, jcenter
        real(wp) dummy

        ! ###################################################################
        ! Get the center
        do j = 1, g%size
            if (g%nodes(j) <= yref .and. &
                g%nodes(j + 1) > yref) then
                jcenter = j
                exit
            end if
        end do

        ! specific potential energy
        if (imode_thermo == THERMO_TYPE_ANELASTIC) then
            ep(:) = (g%nodes - yref)*GRATIO*scaleheightinv
            epbackground(:) = ep(:)
        else
            ep(:) = -(g%nodes - yref)*buoyancy%vector(2)
        end if

        ! hydrstatic pressure
#define p_aux(i)        wrk1d(i,1)
#define r_aux(i)        wrk1d(i,2)

        ! Setting the pressure entry to 1 to get 1/RT
        p_aux(:) = 1.0_wp

        niter = 10

        p(:) = pref             ! initialize iteration
        if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= 0.0_wp) then       ! Get ql, if necessary
            s(:, 3) = 0.0_wp
        end if
        do iter = 1, niter           ! iterate
            if (imode_thermo == THERMO_TYPE_ANELASTIC) then
                pbackground(:) = p_aux(:)
                call THERMO_ANELASTIC_DENSITY(1, g%size, 1, s, r_aux(:))            ! Get r_aux=1/RT
                r_aux(:) = -scaleheightinv*r_aux(:)
            else
                call THERMO_AIRWATER_PH_RE(g%size, s(1, 2), p, s(1, 1), T)
                call THERMO_THERMAL_DENSITY(g%size, s(:, 2), p_aux(:), T, r_aux(:)) ! Get r_aux=1/RT
                r_aux(:) = buoyancy%vector(2)*r_aux(:)
            end if

            p(1) = 0.0_wp
            call OPR_Integral1(1, g, r_aux(:), p, BCS_MIN)

            ! Calculate pressure and normalize s.t. p=pref at y=yref
            p(:) = exp(p(:))
            if (abs(yref - g%nodes(jcenter)) == 0.0_wp) then
                dummy = p(jcenter)
            else
                dummy = p(jcenter) + (p(jcenter + 1) - p(jcenter)) &
                        /(g%nodes(jcenter + 1) - g%nodes(jcenter))*(yref - g%nodes(jcenter))
            end if
            dummy = pref/dummy
            p(:) = dummy*p(:)

            select case (imode_thermo)                  ! Calculate diagnostic variables such as liquid content q_l
            case (THERMO_TYPE_ANELASTIC)
                pbackground(:) = p(:)
                if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= 0.0_wp) then
                    call THERMO_ANELASTIC_PH(1, g%size, 1, s(:, 2), s(:, 1))
                    call THERMO_ANELASTIC_TEMPERATURE(1, g%size, 1, s, T)
                end if

            case (THERMO_TYPE_LINEAR)
                if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
                    call THERMO_AIRWATER_LINEAR(g%size, s, s(:, inb_scal_array))

                end if

            case (THERMO_TYPE_COMPRESSIBLE)
                if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= 0.0_wp) then
                    call THERMO_AIRWATER_PH_RE(g%size, s(1, 2), p, s(1, 1), T)
                end if
            end select

        end do

#undef p_aux
#undef r_aux

        return
    end subroutine Gravity_Hydrostatic_Enthalpy

    !########################################################################
    ! Determine the buoyancy term (density difference rho -rho_0 when it is a function of a scalar
    !########################################################################
    subroutine Gravity_Buoyancy(locProps, nx, ny, nz, s, b, ref)
        type(term_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx, ny, nz, inb_scal_array)
        real(wp), intent(out) :: b(nx, ny, nz)
        real(wp), intent(in) :: ref(ny)         ! reference profile

        ! -----------------------------------------------------------------------
        integer(wi) j, k, is
        real(wp) c0_loc, c1_loc, c2_loc, c3_loc, dummy

        ! #######################################################################
        select case (locProps%type)

        case (EQNS_BOD_HOMOGENEOUS)
            b = locProps%parameters(1)

        case (EQNS_BOD_LINEAR)
            c1_loc = locProps%parameters(1); 
            c2_loc = locProps%parameters(2); 
            c3_loc = locProps%parameters(3)                     ! proportionality factors
            c0_loc = locProps%parameters(inb_scal_array + 1)    ! independent term

            if (locProps%scalar(1) == 1) then
                do k = 1, nz
                    do j = 1, ny
                        dummy = ref(j) - c0_loc
                        b(1:nx, j, k) = c1_loc*s(1:nx, j, k, 1) - dummy
                    end do
                end do

            else if (locProps%scalar(1) == 2) then
                do k = 1, nz
                    do j = 1, ny
                        dummy = ref(j) - c0_loc
                        b(1:nx, j, k) = c1_loc*s(1:nx, j, k, 1) + c2_loc*s(1:nx, j, k, 2) - dummy
                    end do
                end do

            else if (locProps%scalar(1) == 3) then
                do k = 1, nz
                    do j = 1, ny
                        dummy = ref(j) - c0_loc
                        b(1:nx, j, k) = c1_loc*s(1:nx, j, k, 1) + c2_loc*s(1:nx, j, k, 2) + c3_loc*s(1:nx, j, k, 3) - dummy
                    end do
                end do

            else ! general
                do k = 1, nz
                    do j = 1, ny
                        dummy = ref(j)
                        b(1:nx, j, k) = c0_loc - dummy
                    end do
                end do

                do is = 1, locProps%scalar(1)
                    if (abs(locProps%parameters(is)) > small_wp) b(:, :, :) = b(:, :, :) + locProps%parameters(is)*s(:, :, :, is)
                end do

            end if

        case (EQNS_BOD_BILINEAR)
            c0_loc = locProps%parameters(1); 
            c1_loc = locProps%parameters(2); 
            c2_loc = locProps%parameters(3)

            do k = 1, nz
                do j = 1, ny
                    dummy = ref(j)
                    b(1:nx, j, k) = c0_loc*s(1:nx, j, k, 1) + c1_loc*s(1:nx, j, k, 2) + c2_loc*s(1:nx, j, k, 1)*s(1:nx, j, k, 2) - dummy
                end do
            end do

        case (EQNS_BOD_QUADRATIC)
            c0_loc = -locProps%parameters(1)/(locProps%parameters(2)/2.0_wp)**2
            c1_loc = locProps%parameters(2)

            do k = 1, nz
                do j = 1, ny
                    dummy = ref(j)
                    b(1:nx, j, k) = c0_loc*s(1:nx, j, k, 1)*(s(1:nx, j, k, 1) - c1_loc) - dummy
                end do
            end do

        case (EQNS_BOD_NORMALIZEDMEAN)
            c1_loc = locProps%parameters(1)

            do k = 1, nz
                do j = 1, ny
                    dummy = 1.0_wp/bbackground(j)
                    b(1:nx, j, k) = c1_loc*(dummy*s(1:nx, j, k, 1) - 1.0_wp)
                end do
            end do

        case (EQNS_BOD_SUBTRACTMEAN)
            c1_loc = locProps%parameters(1)

            do k = 1, nz
                do j = 1, ny
                    dummy = bbackground(j)
                    b(1:nx, j, k) = c1_loc*(s(1:nx, j, k, 1) - dummy)
                end do
            end do

        case default
            b = 0.0_wp

        end select

        return
    end subroutine Gravity_Buoyancy

    !########################################################################
    !########################################################################
    subroutine Gravity_Buoyancy_Source(locProps, nx, ny, nz, s, gradient, source)
        type(term_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx, ny, nz, inb_scal_array)
        real(wp), intent(inout) :: gradient(nx, ny, nz)         ! gradient magnitude
        real(wp), intent(out) :: source(nx, ny, nz)

        ! -----------------------------------------------------------------------
        real(wp) c0_loc

        ! #######################################################################
        select case (locProps%type)

        case (EQNS_BOD_HOMOGENEOUS)
            source = 0.0_wp

        case (EQNS_BOD_LINEAR)
            source = 0.0_wp

        case (EQNS_BOD_BILINEAR)
            source = 0.0_wp

        case (EQNS_BOD_QUADRATIC)
            c0_loc = -locProps%parameters(1)/(locProps%parameters(2)/2.0_wp)**2; c0_loc = c0_loc*2.0_wp

            source = c0_loc*gradient

        end select

        return
    end subroutine Gravity_Buoyancy_Source

end module Gravity
