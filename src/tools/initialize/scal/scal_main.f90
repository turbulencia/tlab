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
    use Thermodynamics, only: imixture, Thermodynamics_Initialize
    use THERMO_AIRWATER
    use THERMO_ANELASTIC
    use Radiation
    use IO_FIELDS
    use SCAL_LOCAL

    implicit none

! -------------------------------------------------------------------
    integer(wi) is, inb_scal_loc

! ###################################################################
    call TLAB_START()

    call IO_READ_GLOBAL(ifile)
    call Thermodynamics_Initialize(ifile)
    call Radiation_Initialize(ifile)
    call SCAL_READ_LOCAL(ifile)

#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#endif

    if (imode_sim == DNS_MODE_SPATIAL .and. rbg%type == PROFILE_NONE) then
        inb_wrk2d = max(inb_wrk2d, 6)
    end if

    inb_txc = 0
    if (flag_s == PERT_LAYER_BROADBAND) inb_txc = max(inb_txc, 1)
    if (infrared%type /= EQNS_NONE) inb_txc = max(inb_txc, 4)

    call TLAB_ALLOCATE(C_FILE_LOC)

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    call FI_BACKGROUND_INITIALIZE()
    do is = 1, size(IniS)
        if (IniS(is)%relative) IniS(is)%ymean = g(2)%nodes(1) + g(2)%scale*IniS(is)%ymean_rel
    end do

    ! ###################################################################
    itime = 0; rtime = 0.0_wp
    s = 0.0_wp

    ! ###################################################################
    call TLAB_WRITE_ASCII(lfile, 'Initializing scalars.')

    inb_scal_loc = inb_scal
    if (imixture == MIXT_TYPE_AIRWATER) then
        if (damkohler(1) > 0.0_wp .and. flag_mixture == 1) then
            inb_scal_loc = inb_scal - 1
        end if
    end if

    do is = 1, inb_scal_loc
        call SCAL_MEAN(is, s(:, is))

        select case (flag_s)
        case (PERT_LAYER_BROADBAND, PERT_LAYER_DISCRETE)
            call SCAL_FLUCTUATION_VOLUME(is, s(1, is), txc)

        case (PERT_PLANE_BROADBAND, PERT_PLANE_DISCRETE, &
              PERT_DELTA_BROADBAND, PERT_DELTA_DISCRETE, PERT_FLUX_BROADBAND, PERT_FLUX_DISCRETE)
            call SCAL_FLUCTUATION_PLANE(is, s(1, is))

        end select

    end do

    ! ###################################################################
    ! Initial liquid in equilibrium; overwrite previous values
    if (imixture == MIXT_TYPE_AIRWATER) then
        if (damkohler(3) > 0.0_wp .and. flag_mixture == 1) then
            call THERMO_ANELASTIC_PH(imax, jmax, kmax, s(1, 2), s(1, 1))
        end if
    end if

    ! ###################################################################
    ! Initial radiation effect as an accumulation during a certain interval of time
    if (infrared%type /= EQNS_NONE .and. norm_ini_radiation /= 0.0_wp) then

        if (abs(infrared%parameters(1)) > 0.0_wp) then
            infrared%parameters(3) = infrared%parameters(3)/infrared%parameters(1)*norm_ini_radiation
        end if
        infrared%parameters(1) = norm_ini_radiation
        if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= 0.0_wp) then ! Calculate q_l
            call THERMO_ANELASTIC_PH(imax, jmax, kmax, s(1, 2), s(1, 1))
        else if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            call THERMO_AIRWATER_LINEAR(imax*jmax*kmax, s, s(1, inb_scal_array))
        end if
        do is = 1, inb_scal
            if (infrared%active(is)) then
                call Radiation_Infrared_Y(infrared, imax, jmax, kmax, g(2), s, txc(:, 1), txc(:, 2), txc(:, 3), txc(:, 4))
                s(1:isize_field, is) = s(1:isize_field, is) + txc(1:isize_field, 1)
            end if
        end do

    end if

    ! ###################################################################
    call IO_WRITE_FIELDS('scal.ics', IO_SCAL, imax, jmax, kmax, inb_scal, s)

    call TLAB_STOP(0)
end program INISCAL
