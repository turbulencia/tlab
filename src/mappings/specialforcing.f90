#include "dns_const.h"
#include "dns_error.h"

module SpecialForcing
    use TLAB_CONSTANTS, only: wp, wi, pi_wp, efile, MAX_PARS
    use TLAB_TYPES, only: term_dt, grid_dt
    use TLAB_PROCS, only: TLAB_WRITE_ASCII, TLAB_STOP
    implicit none
    private

    type(term_dt) :: forcingProps              ! Forcing parameters

    public :: forcingProps
    public :: SpecialForcing_Initialize
    public :: SpecialForcing_Source
    ! public :: Forcing_Sinusoidal, Forcing_Sinusoidal_NoSlip ! tbd

    integer, parameter :: TYPE_NONE = 0
    integer, parameter :: TYPE_SINUSOIDAL = 1
    integer, parameter :: TYPE_SINUSOIDAL_NOSLIP = 2
    integer, parameter :: TYPE_RAND_MULTIPLICATIVE = 3
    integer, parameter :: TYPE_RAND_ADDIVTIVE = 4
    integer, parameter :: TYPE_WAVEMAKER = 5

contains
    !########################################################################
    !########################################################################
    subroutine SpecialForcing_Initialize(inifile)
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=512) sRes
        integer(wi) idummy
        real(wp) dummy

        !########################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Forcing'

        call TLAB_WRITE_ASCII(bakfile, '#')
        call TLAB_WRITE_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLAB_WRITE_ASCII(bakfile, '#Type=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Geometry=<value>')

        call SCANINICHAR(bakfile, inifile, block, 'Type', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; forcingProps%type = TYPE_NONE
        elseif (trim(adjustl(sRes)) == 'random') then; forcingProps%type = TYPE_RAND_MULTIPLICATIVE; 
        elseif (trim(adjustl(sRes)) == 'sinusoidal') then; forcingProps%type = TYPE_SINUSOIDAL; 
        elseif (trim(adjustl(sRes)) == 'wavemaker') then; forcingProps%type = TYPE_WAVEMAKER; 
        else
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Error in Forcing.Type.')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if

        forcingProps%active(:) = .false.
        if (forcingProps%type /= EQNS_NONE) then
            forcingProps%active(1:3) = .true.       ! default is active in x, y, z momentum equations

            forcingProps%parameters(:) = 0.0_wp
            call SCANINICHAR(bakfile, inifile, block, 'Parameters', '1.0', sRes)
            idummy = MAX_PARS
            call LIST_REAL(sRes, idummy, forcingProps%parameters)

        end if

        forcingProps%vector(:) = 1.0_wp             ! acting equally in the 3 directions
        call SCANINICHAR(bakfile, inifile, block, 'Vector', '0.0,1.0,0.0', sRes)
        idummy = 3
        call LIST_REAL(sRes, idummy, forcingProps%vector)

        select case (forcingProps%type)
        case (TYPE_WAVEMAKER)
            forcingProps%auxiliar(:) = 0.0_wp
            call SCANINICHAR(bakfile, inifile, block, 'Geometry', '1.0,1.0,1.0, 1.0,1.0,1.0', sRes) ! position and region
            idummy = MAX_PARS
            call LIST_REAL(sRes, idummy, forcingProps%auxiliar)

        end select

        ! -------------------------------------------------------------------
        ! backwards compatibility, to be removed
        if (trim(adjustl(sRes)) == 'none') &
            call SCANINIREAL(bakfile, inifile, 'Main', 'TermRandom', '0.0', dummy)
        if (abs(dummy) > 0.0) then
            forcingProps%type = TYPE_RAND_MULTIPLICATIVE
            forcingProps%parameters = dummy
            forcingProps%active(1:3) = .true.
        end if

        return
    end subroutine SpecialForcing_Initialize

!########################################################################
!########################################################################
    subroutine SpecialForcing_Source(locProps, nx, ny, nz, g, time, h, tmp)
        type(term_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz
        type(grid_dt), intent(in) :: g(3)
        real(wp), intent(in) :: time
        real(wp), intent(inout) :: h(nx, ny, nz)
        real(wp), intent(inout) :: tmp(nx, ny, nz)

        ! -----------------------------------------------------------------------
        real(wp) dummy
        integer(wi) i, j, k

        !########################################################################
        select case (locProps%type)

        case (TYPE_RAND_MULTIPLICATIVE)
            dummy = locProps%parameters(1)

            call random_number(tmp)

            tmp = dummy*(tmp*2.0_wp - 1.0_wp)
            tmp = tmp*h

        case (TYPE_SINUSOIDAL)

        case (TYPE_SINUSOIDAL_NOSLIP)

        case (TYPE_WAVEMAKER)
            do j = 1, ny
                do i = 1, nx
                    tmp(i, j, nz) = exp(-0.5_wp*((g(1)%nodes(i) - locProps%auxiliar(1))/locProps%auxiliar(4))**2.0_wp)
                    tmp(i, j, nz) = exp(-0.5_wp*((g(2)%nodes(j) - locProps%auxiliar(2))/locProps%auxiliar(5))**2.0_wp)*tmp(i, j, nz)
                end do
            end do

            if (nz > 1) then
                do k = 1, nz
                    tmp(i, j, k) = exp(0.5_wp*((g(3)%nodes(k) - locProps%auxiliar(3))/locProps%auxiliar(6))**2.0_wp)*tmp(i, j, nz)
                end do
            end if

            tmp = tmp*sin(locProps%parameters(1)*time)

        end select

        return
    end subroutine SpecialForcing_Source

    !########################################################################
    ! Sinusoidal forcing
    !########################################################################
    subroutine Forcing_Sinusoidal(nx, ny, nz, time, visc, u, v, h_u, h_v)
        integer(wi) nx, ny, nz
        real(wp) time, visc
        real(wp), dimension(nx, ny, nz) :: u, v, h_u, h_v

        ! -----------------------------------------------------------------------
        real(wp) omega, sigma, amplitude
        !  integer(wi) i,j

        !########################################################################
        omega = 2.0_wp*pi_wp
        sigma = 2.0_wp*omega*omega*visc

        !  amplitude =-( 1.0_wp + (sigma/omega)**2 )*omega
        !  amplitude = amplitude*sin(omega*time)
        !  DO j = 1,jmax; DO i = 1,imax
        !     u(i,j,1) = SIN(x(i)*omega)*COS(g(2)%nodes(j)*omega)
        !     v(i,j,1) =-COS(x(i)*omega)*SIN(g(2)%nodes(j)*omega)
        !  ENDDO; ENDDO

        amplitude = sin(omega*time)

        h_u = h_u + amplitude*u
        h_v = h_v + amplitude*v

        return
    end subroutine Forcing_Sinusoidal

    !########################################################################
    ! Velocity field with no-slip
    !########################################################################
    subroutine Forcing_Sinusoidal_NoSlip(nx, ny, nz, time, visc, g, h1, h2, tmp1, tmp2, tmp3, tmp4)
        use TLAB_TYPES, only: grid_dt
        use OPR_PARTIAL, only: OPR_PARTIAL_X, OPR_PARTIAL_Y

        integer(wi) nx, ny, nz
        real(wp) time, visc
        type(grid_dt) :: g(:)
        real(wp), dimension(nx*ny*nz) :: h1, h2
        real(wp), dimension(nx*ny*nz) :: tmp1, tmp2, tmp3, tmp4

        ! -----------------------------------------------------------------------
        integer(wi) ij, i, j, bcs(2, 2)

        bcs = 0

        ! #######################################################################
        do j = 1, ny
            do i = 1, nx
                ij = nx*(j - 1) + i
                !     tmp1(ij) = sin(2.0_wp*pi_wp*g(1)%nodes(i))*       sin(C_4_R*pi_wp*g(2)%nodes(j))
                !     tmp2(ij) =-cos(2.0_wp*pi_wp*g(1)%nodes(i))*(1.0_wp-cos(C_4_R*pi_wp*g(2)%nodes(j)))*0.5_wp
                tmp1(ij) = sin(pi_wp*g(1)%nodes(i))*sin(pi_wp*g(1)%nodes(i))*sin(2.0_wp*pi_wp*g(2)%nodes(j))
                tmp2(ij) = -sin(2.0_wp*pi_wp*g(1)%nodes(i))*sin(pi_wp*g(2)%nodes(j))*sin(pi_wp*g(2)%nodes(j))
            end do
        end do

        ! Time terms
        h1 = h1 - tmp1*2.0_wp*pi_wp*sin(2.0_wp*pi_wp*time)
        h2 = h2 - tmp2*2.0_wp*pi_wp*sin(2.0_wp*pi_wp*time)

        ! velocities
        tmp1 = tmp1*cos(2.0_wp*pi_wp*time)
        tmp2 = tmp2*cos(2.0_wp*pi_wp*time)

        ! Diffusion and convection terms in Ox momentum eqn
        call OPR_PARTIAL_Y(OPR_P2_P1, nx, ny, nz, bcs, g(2), tmp1, tmp4, tmp3)
        h1 = h1 - visc*(tmp4) + (tmp3*tmp2)

        call OPR_PARTIAL_X(OPR_P2_P1, nx, ny, nz, bcs, g(1), tmp1, tmp4, tmp3)
        h1 = h1 - visc*(tmp4) + (tmp3*tmp1)

        ! Diffusion and convection terms in Oy momentum eqn
        call OPR_PARTIAL_Y(OPR_P2_P1, nx, ny, nz, bcs, g(2), tmp2, tmp4, tmp3)
        h2 = h2 - visc*(tmp4) + (tmp3*tmp2)

        call OPR_PARTIAL_X(OPR_P2_P1, nx, ny, nz, bcs, g(1), tmp2, tmp4, tmp3)
        h2 = h2 - visc*(tmp4) + (tmp3*tmp1)

        ! #######################################################################
        do j = 1, ny
            do i = 1, nx
                ij = nx*(j - 1) + i
                !     tmp1(ij) = cos(C_4_R*pi_wp*g(1)%nodes(i))*(2.0_wp-cos(C_4_R*pi_wp*g(2)%nodes(j)))/C_8_R &
                !          - 0.5_wp*(sin(2.0_wp*pi_wp*g(2)%nodes(j)))**4
                tmp1(ij) = sin(2.0_wp*pi_wp*g(1)%nodes(i))*sin(2.0_wp*pi_wp*g(2)%nodes(j))
                tmp1(ij) = tmp1(ij)*(cos(2.0_wp*pi_wp*time))**2
            end do
        end do

        ! Pressure gradient
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), tmp1, tmp2)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), tmp1, tmp3)

        h1 = h1 + tmp2
        h2 = h2 + tmp3

        return
    end subroutine Forcing_Sinusoidal_NoSlip

end module SpecialForcing
