#include "dns_error.h"
#include "dns_const.h"

#define C_FILE_LOC "INIFLOW"

program INIFLOW
    use TLAB_CONSTANTS
    use TLAB_VARS
    use TLAB_ARRAYS
    use TLAB_POINTERS, only: e, rho, p, T
    use TLAB_PROCS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_PROCS
#endif
    use Thermodynamics, only: imixture,  Thermodynamics_Initialize
    use THERMO_THERMAL
    use THERMO_CALORIC
    use IO_FIELDS
    use OPR_FOURIER
    use FLOW_LOCAL
    use FLOW_MEAN
    
    implicit none

    !########################################################################
    call TLAB_START()

    call IO_READ_GLOBAL(ifile)
    call Thermodynamics_Initialize(ifile)
    call FLOW_READ_LOCAL(ifile)

#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#endif

    inb_wrk2d = MAX(inb_wrk2d, 3)

    if (flag_u == 0) then
        inb_txc = 2
    else
        inb_txc = 8
    end if

    call TLAB_ALLOCATE(C_FILE_LOC)

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    call FI_BACKGROUND_INITIALIZE()
    if (IniK%relative) IniK%ymean = g(2)%nodes(1) + g(2)%scale*IniK%ymean_rel

    ! Staggering of the pressure grid not implemented here
    if (stagger_on) then
        call TLAB_WRITE_ASCII(wfile, C_FILE_LOC//'. Staggering of the pressure grid not yet implemented.')
        stagger_on = .false. ! turn staggering off for OPR_POISSON_FXZ(...)
    end if
    if (any(PressureFilter%type /= DNS_FILTER_NONE)) then
        call TLAB_WRITE_ASCII(wfile, C_FILE_LOC//'. Pressure and dpdy Filter not implemented here.')
    end if

    if (flag_u /= 0) then ! Initialize Poisson Solver
        if (fourier_on .and. g(1)%periodic .and. g(3)%periodic) then
            call OPR_FOURIER_INITIALIZE()

        else
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. CG routines needed.')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if

    end if

    ! ###################################################################
    itime = 0; rtime = 0.0_wp
    q = 0.0_wp

    ! ###################################################################
    call TLAB_WRITE_ASCII(lfile, 'Initializing velocity.')

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
        call TLAB_WRITE_ASCII(lfile, 'Initializing pressure and density.')

        call PRESSURE_MEAN(p, T, s)
        call DENSITY_MEAN(rho, p, T, s, txc)

        if (flag_u /= 0) then
            call PRESSURE_FLUCTUATION(q(1, 1), q(1, 2), q(1, 3), rho, p, txc(1, 1), &
                                      txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
        end if

        if (imixture > 0) then
            call IO_READ_FIELDS(TRIM(ADJUSTL(tag_scal))//'ics', IO_SCAL, imax, jmax, kmax, inb_scal, 0, s)
        end if

        if (flag_t /= 0) then
            call DENSITY_FLUCTUATION(s, p, rho, txc(1, 1), txc(1, 2))
        end if

        ! Calculate specfic energy. Array s should contain the species fields at this point.
        call THERMO_THERMAL_TEMPERATURE(imax*jmax*kmax, s, p, rho, txc(:, 1))
        call THERMO_CALORIC_ENERGY(imax*jmax*kmax, s, txc(:, 1), e)

    end if

    ! ###################################################################
    call IO_WRITE_FIELDS(TRIM(ADJUSTL(tag_flow))//'ics', IO_FLOW, imax, jmax, kmax, inb_flow, q)

    call TLAB_STOP(0)
end program INIFLOW
