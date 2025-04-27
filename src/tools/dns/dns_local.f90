#include "dns_error.h"
#include "dns_const.h"

module DNS_LOCAL
    use TLab_Constants, only: MAX_VARS, wp, wi, sp
    use TLab_Constants, only: MAX_PATH_LENGTH
#ifdef USE_PSFFT
    use NB3DFFT, only: NB3DFFT_SCHEDLTYPE
#endif
    implicit none
    ! private
    save

    integer :: nitera_first     ! First iteration in current run
    integer :: nitera_last      ! Last iteration in current run
    integer :: nitera_save      ! Iteration step to check-point: save restart files
    integer :: nitera_stats     ! Iteration step to check-point: save statistical data
    integer :: nitera_stats_spa ! Iteration step to accumulate statistics in spatial mode
    integer :: nitera_pln       ! Iteration step to save planes
    integer :: nitera_filter    ! Iteration step for domain filter, if any

    real(wp):: nruntime_sec     ! Maximum runtime of the simulation in seconds
    real(wp):: wall_time        ! Actual elapsed time during the simulation in seconds
    integer :: start_clock      ! Starting time of the simulation on the system

    integer :: nitera_log           ! Iteration step for data logger with simulation information
    character(len=*), parameter :: ofile_base = 'dns.out'    ! data logger filename
    character(len=*), parameter :: vfile_base = 'dns.obs'    ! insitu obs. logger filename
    character(len=MAX_PATH_LENGTH) :: ofile, vfile
    character(len=MAX_PATH_LENGTH) :: logger_path
    real(wp) :: logs_data(20)       ! information (time, time step, cfls, dilatation...)
    real(wp) :: obs_data(20)        ! information (custom variables / insitu measurements ...)
    integer :: dns_obs_log

    integer :: imode_rhs            ! Type of implementation of the RHS of evolution equations
    logical :: remove_divergence    ! Remove residual divergence every time step

    type bounds_dt                      ! control
        sequence
        logical active
        real(wp) min
        real(wp) max
    end type bounds_dt

    type(bounds_dt) bound_p             ! limit pressure in compressible flows
    type(bounds_dt) bound_r             ! limit density in compressible flows
    type(bounds_dt) bound_s(MAX_VARS)    ! limit scalars
    type(bounds_dt) bound_d             ! control dilatation in incompressible/anelastic flows

! Variable viscosity
    logical :: flag_viscosity
    real(wp) :: visc_stop, visc_time, visc_rate

! Tower data (why not in tower module?)
    logical :: use_tower
    integer, dimension(3) :: tower_stride

! NB3DFFT library
#ifdef USE_PSFFT
    type(NB3DFFT_SCHEDLTYPE), save :: nbcsetup
#endif

contains
!########################################################################
! Limit min/max values of fields, if required
!########################################################################
    subroutine DNS_BOUNDS_LIMIT()
        use TLab_Memory, only: inb_scal
        use TLab_Arrays

        ! -------------------------------------------------------------------
        integer(wi) is

        ! ###################################################################
        do is = 1, inb_scal
            if (bound_s(is)%active) then
                s(:, is) = min(max(s(:, is), bound_s(is)%min), bound_s(is)%max)
            end if
        end do

        if (bound_r%active) then
            q(:, 5) = min(max(q(:, 5), bound_r%min), bound_r%max)
        end if

        if (bound_p%active) then
            q(:, 6) = min(max(q(:, 6), bound_p%min), bound_p%max)
        end if

        return
    end subroutine DNS_BOUNDS_LIMIT

!########################################################################
!########################################################################
    subroutine DNS_BOUNDS_CONTROL()
        use TLab_Constants, only: efile, lfile
        use NavierStokes, only: nse_eqns, DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL, DNS_EQNS_ANELASTIC, DNS_EQNS_INCOMPRESSIBLE
        use TLab_WorkFlow, only: stagger_on
        use TLab_Memory, only: imax, jmax, kmax
        use TLab_Arrays
        use TLab_WorkFlow, only: TLab_Write_ASCII
        use Thermo_Anelastic
        use IBM_VARS, only: imode_ibm
#ifdef USE_MPI
        use mpi_f08
        use TLabMPI_VARS, only: ims_offset_i, ims_offset_k
        use TLabMPI_VARS, only: ims_time_min, ims_err
#endif
        use FI_VECTORCALCULUS

        ! -------------------------------------------------------------------
        integer(wi) idummy(3)
        real(wp) dummy
#ifdef USE_MPI
#else
        integer wall_time_loc, int_dummy
#endif
        character*128 line
        character*32 str

        ! Pointers to existing allocated space
        real(wp), dimension(:, :, :), pointer :: loc_max

        ! ###################################################################
        ! Check wall time bounds - maximum runtime
#ifdef USE_MPI
        wall_time = MPI_WTIME() - ims_time_min
        call MPI_BCast(wall_time, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
#else
        ! call ETIME(tdummy, wall_time_loc)
        CALL SYSTEM_CLOCK(wall_time_loc,int_dummy)
        wall_time = real(wall_time_loc - start_clock) / int_dummy
#endif
        ! ###################################################################
        ! Compressible flow
        ! ###################################################################
        select case (nse_eqns)
        case (DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL)

#define p_min_loc logs_data(5)
#define p_max_loc logs_data(6)
#define r_min_loc logs_data(7)
#define r_max_loc logs_data(8)

            ! Check density
            call MINMAX(imax, jmax, kmax, q(:, 5), r_min_loc, r_max_loc)
            if (r_min_loc < bound_r%min .or. r_max_loc > bound_r%max) then
                call TLab_Write_ASCII(efile, 'DNS_CONTROL. Density out of bounds.')
                logs_data(1) = DNS_ERROR_NEGDENS
            end if

            ! Check pressure
            call MINMAX(imax, jmax, kmax, q(:, 6), p_min_loc, p_max_loc)
            if (p_min_loc < bound_p%min .or. p_max_loc > bound_p%max) then
                call TLab_Write_ASCII(efile, 'DNS_CONTROL. Pressure out of bounds.')
                logs_data(1) = DNS_ERROR_NEGPRESS
            end if

        case (DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC)
            if (nse_eqns == DNS_EQNS_ANELASTIC) then
                call Thermo_Anelastic_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, q(1, 1), txc(1, 3))
                call Thermo_Anelastic_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, q(1, 2), txc(1, 4))
                call Thermo_Anelastic_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, q(1, 3), txc(1, 5))
                if (stagger_on) then
                    call FI_INVARIANT_P_STAG(imax, jmax, kmax, txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 1), txc(1, 2), txc(1, 6))
                else
                    call FI_INVARIANT_P(imax, jmax, kmax, txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 1), txc(1, 2))
                end if
            else
                if (stagger_on) then
                    call FI_INVARIANT_P_STAG(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 6))
                else
                    call FI_INVARIANT_P(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2))
                end if
            end if

            if (imode_ibm == 1) then
                if (stagger_on) then
                    call IBM_BCS_FIELD_STAGGER(txc(1, 1)) ! IBM - zeros in solid on pressure mesh
                else
                    call IBM_BCS_FIELD(txc(1, 1))         ! IBM - zeros in solid on velocity mesh
                end if
            end if

#define d_max_loc logs_data(11)
#define d_min_loc logs_data(10)

            call MINMAX(imax, jmax, kmax, txc(1, 1), d_max_loc, d_min_loc)
            d_min_loc = -d_min_loc; d_max_loc = -d_max_loc

            if (max(abs(d_min_loc), abs(d_min_loc)) > bound_d%max) then
                call TLab_Write_ASCII(efile, 'DNS_CONTROL. Dilatation out of bounds.')
                logs_data(1) = DNS_ERROR_DILATATION

                ! Locating the points where the maximum dilatation occurs
                wrk3d = -txc(:, 1)
                loc_max(1:imax, 1:jmax, 1:kmax) => wrk3d(1:imax*jmax*kmax)

                dummy = maxval(wrk3d)
                if (abs(dummy) > bound_d%max) then
                    idummy = maxloc(loc_max)
                    write (str, 1000) dummy; line = 'Maximum dilatation '//trim(adjustl(str))
#ifdef USE_MPI
                    idummy(1) = idummy(1) + ims_offset_i
                    idummy(3) = idummy(3) + ims_offset_k
#endif
                    write (str, *) idummy(1); line = trim(adjustl(line))//' at grid node '//trim(adjustl(str))
                    write (str, *) idummy(2); line = trim(adjustl(line))//':'//trim(adjustl(str))
                    write (str, *) idummy(3); line = trim(adjustl(line))//':'//trim(adjustl(str))//'.'
                    call TLab_Write_ASCII(lfile, line, .true.)
                end if

                dummy = minval(wrk3d)
                if (abs(dummy) > bound_d%max) then
                    idummy = minloc(loc_max)
                    write (str, 1000) dummy; line = 'Minimum dilatation '//trim(adjustl(str))
#ifdef USE_MPI
                    idummy(1) = idummy(1) + ims_offset_i
                    idummy(3) = idummy(3) + ims_offset_k
#endif
                    write (str, *) idummy(1); line = trim(adjustl(line))//' at grid node '//trim(adjustl(str))
                    write (str, *) idummy(2); line = trim(adjustl(line))//':'//trim(adjustl(str))
                    write (str, *) idummy(3); line = trim(adjustl(line))//':'//trim(adjustl(str))//'.'
                    call TLab_Write_ASCII(lfile, line, .true.)
                end if

            end if

        end select

        return

1000    format(G_FORMAT_R)

    end subroutine DNS_BOUNDS_CONTROL
!########################################################################
!########################################################################
    subroutine DNS_OBS_CONTROL()
        use TLab_Memory, only: imax, jmax, kmax, inb_scal
        use TLab_WorkFlow, only: scal_on
        use FDM, only: g
        use FI_VORTICITY_EQN, only: FI_VORTICITY
        use TLab_Arrays
        use Averages
        use Integration, only: Int_Simpson

        integer(wi) :: ip, is

        ! -------------------------------------------------------------------

#define ubulk    obs_data(2)
#define wbulk    obs_data(3)
#define uy1      obs_data(4)
#define wy1      obs_data(5)
#define alpha_1  obs_data(6)
#define alpha_ny obs_data(7)
#define int_ent  obs_data(8)

        ip = 8

        select case (dns_obs_log)

        case (OBS_TYPE_EKMAN)
            ! ubulk, wbulk
            call AVG_IK_V(imax, jmax, kmax, q(1, 1), wrk1d(:, 1), wrk1d(:, 2))
            call AVG_IK_V(imax, jmax, kmax, q(1, 3), wrk1d(:, 3), wrk1d(:, 4))
            ubulk = (1.0_wp/g(2)%nodes(g(2)%size))*Int_Simpson(wrk1d(1:jmax, 1), g(2)%nodes(1:jmax))
            wbulk = (1.0_wp/g(2)%nodes(g(2)%size))*Int_Simpson(wrk1d(1:jmax, 3), g(2)%nodes(1:jmax))

            ! dudy(1), dwdy(1) approximation
            uy1 = wrk1d(2, 1)/g(2)%nodes(2)
            wy1 = wrk1d(2, 3)/g(2)%nodes(2)

            ! turning angles (in degrees)
            alpha_1 = ATAN2D(wy1, uy1)
            alpha_ny = ATAN2D(wrk1d(g(2)%size, 3), wrk1d(g(2)%size, 1))

            ! integrated entstrophy
            call FI_VORTICITY(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3))
            call AVG_IK_V(imax, jmax, kmax, txc(1, 1), wrk1d(:, 1), wrk1d(:, 2))
            int_ent = (1.0_wp/g(2)%nodes(g(2)%size))*Int_Simpson(wrk1d(1:jmax, 1), g(2)%nodes(1:jmax))

            if (scal_on) then
                do is = 1, inb_scal
                    call AVG_IK_V(imax, jmax, kmax, s(1, is), wrk1d(:, 1), wrk1d(:, 2))
                    obs_data(ip + is) = (wrk1d(2, 1) - wrk1d(1, 1))/g(2)%nodes(2)
                end do
            end if

        end select

        return

    end subroutine DNS_OBS_CONTROL

end module DNS_LOCAL
!########################################################################
module DNS_ARRAYS
    use TLab_Constants, only: wp
    implicit none
    save
    private

    real(wp), allocatable, public :: hq(:, :)       ! Right-hand sides Eulerian fields
    real(wp), allocatable, public :: hs(:, :)       ! Right-hand sides Eulerian fields
    real(wp), allocatable, public :: l_hq(:, :)     ! Right-hand sides Lagrangian fields

    target hq, hs, l_hq

end module DNS_ARRAYS
