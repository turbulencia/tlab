#include "dns_error.h"
#include "dns_const.h"

#define C_FILE_LOC "INISCAL"

program INISCAL

    use TLAB_CONSTANTS
    use TLAB_VARS
    use TLAB_ARRAYS
    use TLAB_PROCS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_PROCS
#endif
    use THERMO_VARS, only: imixture
    use THERMO_AIRWATER
    use THERMO_ANELASTIC
    use IO_FIELDS
    use SCAL_LOCAL

    implicit none

! -------------------------------------------------------------------
    integer(wi) is, inb_scal_loc

! ###################################################################
    call TLAB_START()

    call IO_READ_GLOBAL(ifile)
    call THERMO_INITIALIZE()
    call SCAL_READ_LOCAL(ifile)

#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#endif

    isize_wrk3d = isize_field
    if (imode_sim == DNS_MODE_SPATIAL .and. rbg%type == PROFILE_NONE) then
        inb_wrk2d = max(inb_wrk2d, 6)
    end if

    if (flag_s == 1 .or. flag_s == 3 .or. radiation%type /= EQNS_NONE) then; inb_txc = 1
    else; inb_txc = 0
    end if

    call TLAB_ALLOCATE(C_FILE_LOC)

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    call FI_BACKGROUND_INITIALIZE()
    do is = 1, size(Sini)
        if (Sini(is)%relative) Sini(is)%ymean = g(2)%nodes(1) + g(2)%scale*Sini(is)%ymean_rel
    end do

! ###################################################################
    call TLAB_WRITE_ASCII(lfile, 'Initializing scalar fiels.')

    itime = 0; rtime = 0.0_wp
    s = 0.0_wp

    inb_scal_loc = inb_scal
    if (imixture == MIXT_TYPE_AIRWATER) then
        if (damkohler(1) > 0.0_wp .and. flag_mixture == 1) then
            inb_scal_loc = inb_scal - 1
        end if
    end if

    do is = 1, inb_scal_loc
        call SCAL_MEAN(is, s(:, is))

        select case (flag_s)
        case (1, 2, 3)
            call SCAL_FLUCTUATION_VOLUME(is, s(1, is), txc)

        case (4, 5, 6, 7, 8, 9)
            call SCAL_FLUCTUATION_PLANE(is, s(1, is))

        end select

    end do

    if (imixture == MIXT_TYPE_AIRWATER) then ! Initial liquid in equilibrium; overwrite previous values
        if (damkohler(3) > 0.0_wp .and. flag_mixture == 1) then
            call THERMO_ANELASTIC_PH(imax, jmax, kmax, s(1, 2), s(1, 1), epbackground, pbackground)
        end if
    end if

    ! ###################################################################
    if (radiation%type /= EQNS_NONE) then         ! Initial radiation effect as an accumulation during a certain interval of time

        if (ABS(radiation%parameters(1)) > 0.0_wp) then
            radiation%parameters(3) = radiation%parameters(3)/radiation%parameters(1)*norm_ini_radiation
        end if
        radiation%parameters(1) = norm_ini_radiation
        if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= 0.0_wp) then ! Calculate q_l
            call THERMO_ANELASTIC_PH(imax, jmax, kmax, s(1, 2), s(1, 1), epbackground, pbackground)
        else if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            call THERMO_AIRWATER_LINEAR(imax*jmax*kmax, s, s(1, inb_scal_array))
        end if
        do is = 1, inb_scal
            if (radiation%active(is)) then
                call OPR_RADIATION(radiation, imax, jmax, kmax, g(2), s(1, radiation%scalar(is)), txc)
                s(1:isize_field, is) = s(1:isize_field, is) + txc(1:isize_field, 1)
            end if
        end do

    end if

! ###################################################################
    call IO_WRITE_FIELDS('scal.ics', IO_SCAL, imax, jmax, kmax, inb_scal, s)

    call TLAB_STOP(0)
end program INISCAL
