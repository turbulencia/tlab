#include "dns_error.h"
#include "dns_const.h"

#define C_FILE_LOC "INISCAL"

program INISCAL

    use TLab_Constants
    use TLAB_VARS
    use TLab_Arrays
    use TLab_WorkFlow
    use TLab_Memory, only: TLab_Initialize_Memory
#ifdef USE_MPI
    use MPI
    use TLabMPI_PROCS
#endif
    use FDM, only: g,  FDM_Initialize
    use Thermodynamics, only: imixture, Thermodynamics_Initialize_Parameters
    use THERMO_AIRWATER
    use THERMO_ANELASTIC
    use Radiation
    use IO_FIELDS
    use SCAL_LOCAL

    implicit none

! -------------------------------------------------------------------
    integer(wi) is, inb_scal_loc

! ###################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize()
#endif
    call NavierStokes_Initialize_Parameters(ifile)
    call Thermodynamics_Initialize_Parameters(ifile)

    call SCAL_READ_LOCAL(ifile)

    if (imode_sim == DNS_MODE_SPATIAL) then
        inb_wrk2d = max(inb_wrk2d, 6)
    end if

    inb_txc = 0
    if (flag_s == PERT_LAYER_BROADBAND) inb_txc = max(inb_txc, 1)
    if (norm_ini_radiation /= 0.0_wp) inb_txc = max(inb_txc, 4)

    call TLab_Initialize_Memory(C_FILE_LOC)

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, wrk1d(:,1), wrk1d(:,2), wrk1d(:,3))
    call FDM_INITIALIZE(x, g(1), wrk1d(:,1), wrk1d(:,4))
    call FDM_INITIALIZE(y, g(2), wrk1d(:,2), wrk1d(:,4))
    call FDM_INITIALIZE(z, g(3), wrk1d(:,3), wrk1d(:,4))

    call Radiation_Initialize(ifile)

    call TLab_Initialize_Background()
    do is = 1, size(IniS)
        if (IniS(is)%relative) IniS(is)%ymean = g(2)%nodes(1) + g(2)%scale*IniS(is)%ymean_rel
    end do

    ! ###################################################################
    itime = 0; rtime = 0.0_wp
    s = 0.0_wp

    ! ###################################################################
    call TLab_Write_ASCII(lfile, 'Initializing scalars.')

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
    if (infraredProps%type /= EQNS_NONE .and. norm_ini_radiation /= 0.0_wp) then
        norm_ini_radiation = norm_ini_radiation/infraredProps%auxiliar(1)
        infraredProps%auxiliar(:) = infraredProps%auxiliar(:)*norm_ini_radiation
        if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= 0.0_wp) then ! Calculate q_l
            call THERMO_ANELASTIC_PH(imax, jmax, kmax, s(1, 2), s(1, 1))
        else if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            call THERMO_AIRWATER_LINEAR(imax*jmax*kmax, s, s(1, inb_scal_array))
        end if
        do is = 1, inb_scal
            if (infraredProps%active(is)) then
                call Radiation_Infrared_Y(infraredProps, imax, jmax, kmax, g(2), s, txc(:, 1), txc(:, 2), txc(:, 3), txc(:, 4))
                s(1:isize_field, is) = s(1:isize_field, is) + txc(1:isize_field, 1)
            end if
        end do

    end if

    ! ###################################################################
    call IO_WRITE_FIELDS('scal.ics', IO_SCAL, imax, jmax, kmax, inb_scal, s)

    call TLab_Stop(0)
end program INISCAL
