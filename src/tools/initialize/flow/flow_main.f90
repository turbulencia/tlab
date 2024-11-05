#include "dns_error.h"
#include "dns_const.h"

#define C_FILE_LOC "INIFLOW"

program INIFLOW
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: ifile, gfile, lfile, efile, wfile, tag_flow, tag_scal
    use TLAB_VARS
    use TLab_Arrays
    use TLab_Pointers, only: e, rho, p, T
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Memory, only: TLab_Initialize_Memory
#ifdef USE_MPI
    use MPI
    use TLabMPI_PROCS
#endif
    use FDM, only: g, x, y, z, FDM_Initialize
    use Thermodynamics, only: imixture, Thermodynamics_Initialize_Parameters
    use THERMO_THERMAL
    use THERMO_CALORIC
    use IO_FIELDS
    use OPR_Elliptic
    use OPR_FOURIER
    use FLOW_LOCAL
    use FLOW_MEAN

    implicit none

    !########################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize()
#endif
    call NavierStokes_Initialize_Parameters(ifile)
    call Thermodynamics_Initialize_Parameters(ifile)
    call FLOW_READ_LOCAL(ifile)

    inb_wrk2d = max(inb_wrk2d, 3)

    if (flag_u == 0) then
        inb_txc = 2
    else
        inb_txc = 8
    end if

    call TLab_Initialize_Memory(C_FILE_LOC)

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, wrk1d(:,1), wrk1d(:,2), wrk1d(:,3))
    call FDM_INITIALIZE(x, g(1), wrk1d(:,1), wrk1d(:,4))
    call FDM_INITIALIZE(y, g(2), wrk1d(:,2), wrk1d(:,4))
    call FDM_INITIALIZE(z, g(3), wrk1d(:,3), wrk1d(:,4))

    call TLab_Initialize_Background()
    if (IniK%relative) IniK%ymean = g(2)%nodes(1) + g(2)%scale*IniK%ymean_rel

    ! Staggering of the pressure grid not implemented here
    if (stagger_on) then
        call TLab_Write_ASCII(wfile, C_FILE_LOC//'. Staggering of the pressure grid not yet implemented.')
        stagger_on = .false. ! turn staggering off for OPR_Poisson_FourierXZ_Factorize(...)
    end if
    if (any(PressureFilter%type /= DNS_FILTER_NONE)) then
        call TLab_Write_ASCII(wfile, C_FILE_LOC//'. Pressure and dpdy Filter not implemented here.')
    end if

    if (flag_u /= 0) then ! Initialize Poisson Solver
        call OPR_Elliptic_Initialize(ifile)

        if (fourier_on .and. g(1)%periodic .and. g(3)%periodic) then
            call OPR_FOURIER_INITIALIZE()

        else
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. CG routines needed.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

    end if

    ! ###################################################################
    itime = 0; rtime = 0.0_wp
    q = 0.0_wp

    ! ###################################################################
    call TLab_Write_ASCII(lfile, 'Initializing velocity.')

    call VELOCITY_MEAN(q(1, 1), q(1, 2), q(1, 3))

    select case (flag_u)
    case (PERT_DISCRETE)
        call VELOCITY_DISCRETE(txc(1, 1), txc(1, 2), txc(1, 3))
        q(1:isize_field, 1:3) = q(1:isize_field, 1:3) + txc(1:isize_field, 1:3)

    case (PERT_BROADBAND, PERT_BROADBAND_POTENTIAL, PERT_BROADBAND_VORTICITY)
        call VELOCITY_BROADBAND(txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7), txc(1, 8))
        q(1:isize_field, 1:3) = q(1:isize_field, 1:3) + txc(1:isize_field, 1:3)

    end select

    ! ###################################################################
    ! Compressible formulation
    if (imode_eqns == DNS_EQNS_TOTAL .or. imode_eqns == DNS_EQNS_INTERNAL) then
        call TLab_Write_ASCII(lfile, 'Initializing pressure and density.')

        call PRESSURE_MEAN(p, T, s)
        call DENSITY_MEAN(rho, p, T, s, txc)

        if (flag_u /= 0) then
            call PRESSURE_FLUCTUATION(q(1, 1), q(1, 2), q(1, 3), rho, p, txc(1, 1), &
                                      txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
        end if

        if (imixture > 0) then
            call IO_READ_FIELDS(trim(adjustl(tag_scal))//'ics', IO_SCAL, imax, jmax, kmax, inb_scal, 0, s)
        end if

        if (flag_t /= 0) then
            call DENSITY_FLUCTUATION(s, p, rho, txc(1, 1), txc(1, 2))
        end if

        ! Calculate specfic energy. Array s should contain the species fields at this point.
        call THERMO_THERMAL_TEMPERATURE(imax*jmax*kmax, s, p, rho, txc(:, 1))
        call THERMO_CALORIC_ENERGY(imax*jmax*kmax, s, txc(:, 1), e)

    end if

    ! ###################################################################
    call IO_WRITE_FIELDS(trim(adjustl(tag_flow))//'ics', IO_FLOW, imax, jmax, kmax, inb_flow, q)

    call TLab_Stop(0)
end program INIFLOW
