#include "dns_error.h"
#include "dns_const.h"

#define C_FILE_LOC "INIFLOW"

program INIFLOW
    use TLAB_CONSTANTS
    use TLAB_VARS
    use TLAB_ARRAYS
    use TLAB_PROCS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_PROCS
#endif
    use THERMO_VARS, only: imixture
    use FLOW_LOCAL
    use IO_FIELDS
#ifdef USE_CGLOC
    use CG_GLOBAL, only: cg_unif, cg_ord
#endif

    implicit none

    ! -------------------------------------------------------------------
    ! Additional local arrays
#ifdef USE_CGLOC
    real(wp), dimension(:), allocatable, save :: ci, cj, ck, ipos, jpos, kpos
#endif

    !########################################################################
    call TLAB_START()

    call IO_READ_GLOBAL(ifile)
    call FLOW_READ_LOCAL(ifile)

#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#endif

    inb_wrk2d = MAX(inb_wrk2d, 3)
    isize_wrk3d = isize_txc_field

    if (flag_u == 0) then; inb_txc = 2
    else; inb_txc = 8
    end if

    call TLAB_ALLOCATE(C_FILE_LOC)

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    call FI_BACKGROUND_INITIALIZE(wrk1d)
    if (Kini%relative) Kini%ymean = g(2)%nodes(1) + g(2)%scale*Kini%ymean_rel

    ! Staggering of the pressure grid not implemented here
    if (istagger == 1 .or. ivfilter == 1) then
        call TLAB_WRITE_ASCII(wfile, C_FILE_LOC//'. Staggering of the pressure grid not yet implemented.')
        istagger = 0; ivfilter = 0 ! turn staggering off for OPR_POISSON_FXZ(...)
    end if

    if (flag_u /= 0) then ! Initialize Poisson Solver
        if (ifourier == 1 .and. g(1)%periodic .and. g(3)%periodic) then
            call OPR_FOURIER_INITIALIZE(txc, wrk1d, wrk2d, wrk3d)

        else
#ifdef USE_CGLOC
            allocate (ci(isize_field*2))
            allocate (cj(isize_field*2))
            allocate (ck(isize_field*2))
            allocate (ipos(imax*4))
            allocate (jpos(jmax*4))
            allocate (kpos(kmax*4))

            if (.not. g(1)%uniform .or. .not. g(2)%uniform) then
                call TLAB_WRITE_ASCII(lfile, 'Initializing conjugate gradient, non-uniform grid, second-order.')
                cg_unif = 1; cg_ord = 2
                ! to be rewritten in terms of grid derived type
                ! CALL CGBC2(cg_unif, imode_fdm, imax,jmax,kmax,g(3)%size, &
                !      i1bc,j1bc,k1bc, scalex,scaley,scalez, dx,dy,dz, ipos,jpos,kpos,ci,cj,ck, wrk2d)
            else
                call TLAB_WRITE_ASCII(lfile, 'Initializing conjugate gradient, uniform grid, fourth-order.')
                cg_unif = 0; cg_ord = 4
                ! CALL CGBC4(cg_unif, imax,jmax,kmax,g(3)%size, &
                !      i1bc,j1bc,k1bc, scalex,scaley,scalez, dx,dy,dz, ipos,jpos,kpos,ci,cj,ck, wrk2d)
            end if
#else
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. CG routines needed.')
            call TLAB_STOP(DNS_ERROR_OPTION)
#endif
        end if

    end if

    ! ###################################################################
    itime = 0; rtime = 0.0_wp
    q = 0.0_wp

    ! ###################################################################
    call TLAB_WRITE_ASCII(lfile, 'Initializing velocity.')

    call VELOCITY_MEAN(q(1, 1), q(1, 2), q(1, 3), wrk1d)

    select case (flag_u)
    case (1)
        call VELOCITY_DISCRETE(txc(1, 1), txc(1, 2), txc(1, 3), wrk1d, wrk2d)
        q(1:isize_field, 1:3) = q(1:isize_field, 1:3) + txc(1:isize_field, 1:3)

    case (2, 3, 4)
        call VELOCITY_BROADBAND(txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7), txc(1, 8), &
                                wrk1d, wrk2d, wrk3d)
        q(1:isize_field, 1:3) = q(1:isize_field, 1:3) + txc(1:isize_field, 1:3)

    end select

    ! ###################################################################
    if (imode_eqns == DNS_EQNS_TOTAL .or. imode_eqns == DNS_EQNS_INTERNAL) then
        call TLAB_WRITE_ASCII(lfile, 'Initializing pressure and density.')

        ! for readibility: energy, density, pressure and temperature
#define e_loc      q(:, 4)
#define r_loc      q(:, 5)
#define p_loc      q(:, 6)
#define T_loc      q(:, 7)
        call PRESSURE_MEAN(p_loc, T_loc, s, wrk1d)
        call DENSITY_MEAN(r_loc, p_loc, T_loc, s, txc, wrk1d, wrk2d, wrk3d)

        if (flag_u /= 0) then
            call PRESSURE_FLUCTUATION(q(1, 1), q(1, 2), q(1, 3), r_loc, p_loc, txc(1, 1), &
                                      txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), wrk1d, wrk2d, wrk3d)
        end if

        if (imixture > 0) then
            call IO_READ_FIELDS(TRIM(ADJUSTL(tag_scal))//'ics', IO_SCAL, imax, jmax, kmax, inb_scal, 0, s, wrk3d)
        end if

        if (flag_t == 4 .or. flag_t == 5) then
            call DENSITY_FLUCTUATION(flag_t, s, p_loc, r_loc, txc(1, 1), txc(1, 2), wrk2d)
        end if

        ! Calculate specfic energy. Array s should contain the species fields at this point.
        call THERMO_THERMAL_TEMPERATURE(imax, jmax, kmax, s, p_loc, r_loc, txc(1, 1))
        call THERMO_CALORIC_ENERGY(imax, jmax, kmax, s, txc(1, 1), e_loc)

    end if

    ! ###################################################################
    call IO_WRITE_FIELDS(TRIM(ADJUSTL(tag_flow))//'ics', IO_FLOW, imax, jmax, kmax, inb_flow, q, wrk3d)

    call TLAB_STOP(0)
end program INIFLOW
