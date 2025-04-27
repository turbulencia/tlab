#include "dns_const.h"
#include "dns_error.h"

module SpecialForcing
    use TLab_Constants, only: wp, wi, pi_wp, efile, lfile, MAX_VARS, MAX_PARS
    use FDM, only: fdm_dt
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Memory, only: TLab_Allocate_Real
    use TLab_Arrays, only: wrk1d
    implicit none
    private

    type term_dt
        sequence
        integer type
        integer scalar(MAX_VARS)                ! fields defining this term
        logical active(MAX_VARS), lpadding(3)   ! fields affected by this term
        real(wp) parameters(MAX_PARS)
        real(wp) auxiliar(MAX_PARS)
        real(wp) vector(3)
    end type term_dt
    type(term_dt), public, protected :: forcingProps              ! Forcing parameters

    public :: SpecialForcing_Initialize
    public :: SpecialForcing_Source
    ! public :: Forcing_Sinusoidal, Forcing_Sinusoidal_NoSlip ! tbd

    integer, parameter :: TYPE_NONE = 0
    integer, parameter :: TYPE_HOMOGENEOUS = 1
    integer, parameter :: TYPE_SINUSOIDAL = 2
    integer, parameter :: TYPE_SINUSOIDAL_NOSLIP = 3
    integer, parameter :: TYPE_RAND_MULTIPLICATIVE = 4
    integer, parameter :: TYPE_RAND_ADDIVTIVE = 5
    integer, parameter :: TYPE_WAVEMAKER = 6

    integer, parameter :: nwaves_max = 3                ! maximum number of waves in wavemaker
    integer :: nwaves                                   ! number of waves in wavemaker
    real(wp) amplitude(2, nwaves_max)                   ! wave amplitudes in x, y
    real(wp) wavenumber(2, nwaves_max)                  ! wavenumbers in x, y
    real(wp) frequency(nwaves_max)                      ! wave frequencies
    real(wp) envelope(4)                                ! (x,y,z) position, size

    real(wp), allocatable, target :: tmp_envelope(:, :) ! arrays for this routine
    real(wp), allocatable, target :: tmp_phase(:, :)

contains
    !########################################################################
    !########################################################################
    subroutine SpecialForcing_Initialize(inifile)
        use TLab_Memory, only: imax, jmax, kmax
#ifdef USE_MPI
        use TLabMPI_VARS, only: ims_offset_i, ims_offset_k
#endif
        use FDM, only: g
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=512) sRes
        integer(wi) idummy, i, j, k, iwave, idsp, kdsp
        real(wp) :: dummy(MAX_PARS)

        real(wp), pointer :: p_envelope(:, :, :) => null(), p_phase(:, :, :) => null()

        !########################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'SpecialForcing'

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<value>')
        call TLab_Write_ASCII(bakfile, '#Parameters=<values>')
        call TLab_Write_ASCII(bakfile, '#Wave#=<amplitude,wavenumber,angle,frequency>')
        call TLab_Write_ASCII(bakfile, '#Envelope=<x,y,z,size>')

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; forcingProps%type = TYPE_NONE
        elseif (trim(adjustl(sRes)) == 'homogeneous') then; forcingProps%type = TYPE_HOMOGENEOUS; 
        elseif (trim(adjustl(sRes)) == 'random') then; forcingProps%type = TYPE_RAND_MULTIPLICATIVE; 
        elseif (trim(adjustl(sRes)) == 'sinusoidal') then; forcingProps%type = TYPE_SINUSOIDAL; 
        elseif (trim(adjustl(sRes)) == 'wavemaker') then; forcingProps%type = TYPE_WAVEMAKER; 
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Error in SpecialForcing.Type.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        forcingProps%vector(:) = 0.0_wp; forcingProps%active(:) = .false.
        if (forcingProps%type /= TYPE_NONE) then
            call ScanFile_Char(bakfile, inifile, block, 'Vector', '1.0, 0.0, 0.0', sRes)
            idummy = 3
            call LIST_REAL(sRes, idummy, forcingProps%vector)

            if (abs(forcingProps%vector(1)) > 0.0_wp) then; forcingProps%active(1) = .true.; call TLab_Write_ASCII(lfile, 'Forcing along Ox.'); end if
            if (abs(forcingProps%vector(2)) > 0.0_wp) then; forcingProps%active(2) = .true.; call TLab_Write_ASCII(lfile, 'Forcing along Oy.'); end if
            if (abs(forcingProps%vector(3)) > 0.0_wp) then; forcingProps%active(3) = .true.; call TLab_Write_ASCII(lfile, 'Forcing along Oz.'); end if

            forcingProps%parameters(:) = 0.0_wp
            call ScanFile_Char(bakfile, inifile, block, 'Parameters', '1.0, 1.0, 0.0', sRes)
            idummy = MAX_PARS
            call LIST_REAL(sRes, idummy, forcingProps%parameters)

            select case (forcingProps%type)
            case (TYPE_HOMOGENEOUS)

            case (TYPE_WAVEMAKER)
                forcingProps%active(1:3) = .true.       ! default is active in x, y, z momentum equations

                do nwaves = 1, nwaves_max
                    write (sRes, *) nwaves
                    call ScanFile_Char(bakfile, inifile, block, 'Wave'//trim(adjustl(sRes)), 'void', sRes)
                    if (trim(adjustl(sRes)) /= 'void') then
                        idummy = 4
                        call LIST_REAL(sRes, idummy, dummy)
                        if (idummy /= 4) then
                            call TLab_Write_ASCII(efile, __FILE__//'. Error in '//trim(adjustl(block))//'.Wave.')
                            call TLab_Stop(DNS_ERROR_OPTION)
                        end if
                        dummy(3) = dummy(3)*pi_wp/180._wp                   ! from degree to radians
                        wavenumber(1, nwaves) = dummy(2)*cos(dummy(3))      ! x-wavenumber
                        wavenumber(2, nwaves) = dummy(2)*sin(dummy(3))      ! y-wavenumber
                        amplitude(1, nwaves) = dummy(1)*sin(dummy(3))       ! in x equation
                        amplitude(2, nwaves) = -dummy(1)*cos(dummy(3))      ! in y equation
                        frequency(nwaves) = dummy(4)
                    else
                        exit
                    end if
                end do
                nwaves = nwaves - 1                                         ! correct for the increment in the loop

                envelope(:) = 0.0_wp
                call ScanFile_Char(bakfile, inifile, block, 'Envelope', '1.0,1.0,1.0, 1.0', sRes) ! position and size
                idummy = MAX_PARS
                call LIST_REAL(sRes, idummy, envelope)
                envelope(4) = abs(envelope(4))                              ! make sure the size parameter is positive

                forcingProps%active(3) = .false.                            ! only active in x and y

            end select

        end if

        !########################################################################
        select case (forcingProps%type)
        case (TYPE_WAVEMAKER)
            call TLab_Allocate_Real(__FILE__, tmp_envelope, [imax*jmax, kmax], 'tmp-wave-envelope')
            call TLab_Allocate_Real(__FILE__, tmp_phase, [imax*jmax, nwaves], 'tmp-wave-phase')

#ifdef USE_MPI
            idsp = ims_offset_i; kdsp = ims_offset_k
#else
            idsp = 0; kdsp = 0
#endif
            p_envelope(1:imax, 1:jmax, 1:kmax) => tmp_envelope(1:imax*jmax*kmax, 1)
            p_phase(1:imax, 1:jmax, 1:nwaves) => tmp_phase(1:imax*jmax*nwaves, 1)

            dummy(1) = 0.5_wp/envelope(4)**2.0_wp
            wrk1d(1:imax, 1) = g(1)%nodes(idsp + 1:idsp + imax) - envelope(1)
            wrk1d(1:jmax, 2) = g(2)%nodes(1:jmax) - envelope(2)
            wrk1d(1:kmax, 3) = g(3)%nodes(kdsp + 1:kdsp + kmax) - envelope(3)
            do k = 1, kmax
                do j = 1, jmax
                    do i = 1, imax
                        p_envelope(i, j, k) = wrk1d(i, 1)**2.0_wp + wrk1d(j, 2)**2.0_wp + wrk1d(k, 3)**2.0_wp
                        p_envelope(i, j, k) = exp(-dummy(1)*p_envelope(i, j, k))         ! exp of an array can cause memory problems
                        do iwave = 1, nwaves
                            p_phase(i, j, iwave) = wrk1d(i, 1)*wavenumber(1, iwave) + wrk1d(j, 2)*wavenumber(2, iwave)
                        end do
                    end do
                end do
            end do

            nullify (p_envelope, p_phase)

        end select

        ! -------------------------------------------------------------------
        ! Check with previous version; to be removed
        call ScanFile_Char(bakfile, inifile, 'Main', 'TermRandom', 'void', sRes)
        if (trim(adjustl(sRes)) /= 'void') then
            call TLab_Write_ASCII(efile, __FILE__//'. Update TermRandom to [SpecialForcing].')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        return
    end subroutine SpecialForcing_Initialize

!########################################################################
!########################################################################
    subroutine SpecialForcing_Source(locProps, nx, ny, nz, iq, time, q, h, tmp)
        type(term_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz, iq
        real(wp), intent(in) :: time
        real(wp), intent(in) :: q(nx*ny, nz)
        real(wp), intent(inout) :: h(nx*ny, nz)
        real(wp), intent(inout) :: tmp(nx*ny, nz)

        ! -----------------------------------------------------------------------
        integer(wi) k, iwave

        !########################################################################
        select case (locProps%type)

        case (TYPE_HOMOGENEOUS)
            tmp(:, :) = locProps%parameters(1)

        case (TYPE_RAND_MULTIPLICATIVE)
            call random_number(tmp)

            tmp(:, :) = (tmp(:, :)*2.0_wp - 1.0_wp)*locProps%parameters(1)
            tmp(:, :) = tmp(:, :)*h(:, :)

        case (TYPE_SINUSOIDAL)

        case (TYPE_SINUSOIDAL_NOSLIP)

        case (TYPE_WAVEMAKER)
            tmp(:, :) = 0.0_wp
            do k = 1, nz
                do iwave = 1, nwaves
                    tmp(:, k) = tmp(:, k) + sin(tmp_phase(:, iwave) - frequency(iwave)*time)*amplitude(iq, iwave)
                end do
                tmp(:, k) = (tmp(:, k) - q(:, k))*tmp_envelope(:, k)*locProps%parameters(1)
            end do

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
        use FDM, only: fdm_dt
        use OPR_Partial, only: OPR_Partial_X, OPR_Partial_Y, OPR_P1, OPR_P2_P1

        integer(wi) nx, ny, nz
        real(wp) time, visc
        type(fdm_dt) :: g(:)
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
        call OPR_Partial_Y(OPR_P2_P1, nx, ny, nz, bcs, g(2), tmp1, tmp4, tmp3)
        h1 = h1 - visc*(tmp4) + (tmp3*tmp2)

        call OPR_Partial_X(OPR_P2_P1, nx, ny, nz, bcs, g(1), tmp1, tmp4, tmp3)
        h1 = h1 - visc*(tmp4) + (tmp3*tmp1)

        ! Diffusion and convection terms in Oy momentum eqn
        call OPR_Partial_Y(OPR_P2_P1, nx, ny, nz, bcs, g(2), tmp2, tmp4, tmp3)
        h2 = h2 - visc*(tmp4) + (tmp3*tmp2)

        call OPR_Partial_X(OPR_P2_P1, nx, ny, nz, bcs, g(1), tmp2, tmp4, tmp3)
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
        call OPR_Partial_X(OPR_P1, nx, ny, nz, bcs, g(1), tmp1, tmp2)
        call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs, g(2), tmp1, tmp3)

        h1 = h1 + tmp2
        h2 = h2 + tmp3

        return
    end subroutine Forcing_Sinusoidal_NoSlip

end module SpecialForcing
