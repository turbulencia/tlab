#include "dns_const.h"
#include "dns_error.h"
#include "avgij_map.h"

module STATISTICS
    use TLAB_CONSTANTS, only: MAX_AVG_TEMPORAL, wp, wi, small_wp
    implicit none
    save
    private

    real(wp), allocatable, public :: mean(:, :)
    real(wp), allocatable, public :: mean_flow(:, :, :)     ! These 2 are for spatial case
    real(wp), allocatable, public :: mean_scal(:, :, :, :)

    logical, public :: stats_averages, stats_pdfs, stats_intermittency, stats_buoyancy

    public :: STATISTICS_INITIALIZE, STATISTICS_TEMPORAL, STATISTICS_SPATIAL

contains

    ! ###################################################################
    ! ###################################################################
    subroutine STATISTICS_INITIALIZE()

        use TLAB_VARS, only: imode_sim, jmax, inb_scal, nstatavg

        if (imode_sim == DNS_MODE_TEMPORAL) then
            allocate (mean(jmax, MAX_AVG_TEMPORAL))

        else if (imode_sim == DNS_MODE_SPATIAL) then
            allocate (mean_flow(nstatavg, jmax, MA_MOMENTUM_SIZE))
            allocate (mean_scal(nstatavg, jmax, MS_SCALAR_SIZE, inb_scal))

        end if

        return
    end subroutine STATISTICS_INITIALIZE

    !########################################################################
    !########################################################################
    subroutine STATISTICS_TEMPORAL()

#ifdef TRACE_ON
        use TLAB_CONSTANTS, only: tfile
        use TLAB_PROCS, only: TLAB_WRITE_ASCII
#endif
        use TLAB_TYPES, only: pointers_dt
        use TLAB_VARS, only: g
        use TLAB_VARS, only: imax, jmax, kmax, isize_field, inb_scal_array
        use TLAB_VARS, only: buoyancy, imode_eqns, scal_on
        use TLAB_VARS, only: froude
        use TLAB_VARS, only: itime, rtime
        use TLAB_VARS, only: schmidt
        use TLAB_ARRAYS
        use THERMO_ANELASTIC
        use DNS_ARRAYS
        use Thermodynamics, only: imixture
        use PARTICLE_VARS
        use PARTICLE_ARRAYS
        use FI_SOURCES, only: FI_BUOYANCY
        use FI_VORTICITY_EQN

        ! -------------------------------------------------------------------
        real(wp) dummy, amin(16), amax(16)
        integer is, idummy, nbins, ibc(16), nfield
        integer(wi) ij
        type(pointers_dt) vars(16)
        character*32 fname, gatename(1)
        character*64 str
        integer(1) igate
        integer(1), allocatable, save :: gate(:)

        ! ###################################################################
#ifdef TRACE_ON
        call TLAB_WRITE_ASCII(tfile, 'ENTERING STATS_TEMPORAL_LAYER')
#endif

        stats_buoyancy = .false.  ! default

        ! in case we need the buoyancy statistics
        if (buoyancy%type == EQNS_BOD_QUADRATIC .or. &
            buoyancy%type == EQNS_BOD_BILINEAR .or. &
            imixture == MIXT_TYPE_AIRWATER .or. &
            imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            stats_buoyancy = .true.
        end if

        ! Calculate pressure
        if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
            call FI_PRESSURE_BOUSSINESQ(q, s, txc(1, 3), txc(1, 1), txc(1, 2), txc(1, 4))
        end if

        ! ###################################################################
        ! Intermittency
        ! ###################################################################
        if (stats_intermittency) then
            allocate (gate(isize_field))
            call FI_VORTICITY(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 4))

            ! calculate vorticity gate based on 0.1% threshold
            call MINMAX(imax, jmax, kmax, txc(1, 1), amin(1), amax(1))
            amin(1) = 1.0e-3_wp*1.0e-3_wp*amax(1)
            do ij = 1, isize_field
                if (txc(ij, 1) > amin(1)) then; gate(ij) = 1  ! gate array
                else; gate(ij) = 0
                end if
            end do
            nfield = 1; gatename(1) = 'Vorticity'

            write (fname, *) itime; fname = 'int'//trim(adjustl(fname))
            call INTER_N_XZ(fname, itime, rtime, imax, jmax, kmax, nfield, gatename, gate, g(2)%nodes, mean)

            deallocate (gate)
        end if

        ! ###################################################################
        ! Unconditional plane PDFs
        ! ###################################################################
        if (stats_pdfs) then
            nfield = 0
            nfield = nfield + 1; vars(nfield)%field => q(:, 1); vars(nfield)%tag = 'u'
            nfield = nfield + 1; vars(nfield)%field => q(:, 2); vars(nfield)%tag = 'v'
            nfield = nfield + 1; vars(nfield)%field => q(:, 3); vars(nfield)%tag = 'w'
            if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
                nfield = nfield + 1; vars(nfield)%field => txc(:, 3); vars(nfield)%tag = 'p'
            else
                nfield = nfield + 1; vars(nfield)%field => q(:, 6); vars(nfield)%tag = 'p'
                nfield = nfield + 1; vars(nfield)%field => q(:, 5); vars(nfield)%tag = 'r'
                nfield = nfield + 1; vars(nfield)%field => q(:, 7); vars(nfield)%tag = 't'
            end if

            do is = 1, inb_scal_array
                nfield = nfield + 1; vars(nfield)%field => s(:, is); vars(nfield)%tag = 's'
                write (str, *) is; vars(nfield)%tag = trim(adjustl(vars(nfield)%tag))//trim(adjustl(str))
            end do

            ibc(1:nfield) = 2 ! BCs in the calculation of the PDFs
            igate = 0         ! no intermittency partition

            nbins = 32
            write (fname, *) itime; fname = 'pdf'//trim(adjustl(fname))
            call PDF1V_N(fname, rtime, imax, jmax, kmax, &
                         nfield, nbins, ibc, amin, amax, vars, igate, wrk3d, g(2)%nodes, txc)

        end if

        ! ###################################################################
        ! Plane averages
        ! ###################################################################
        if (stats_averages) then
            if (scal_on) then
                do is = 1, inb_scal_array          ! All, prognostic and diagnostic fields in array s
                    hq(1:isize_field, 3) = txc(1:isize_field, 3) ! Pass the pressure
                    call AVG_SCAL_XZ(is, q, s, s(1, is), &
                                     txc(1, 1), txc(1, 2), txc(1, 4), txc(1, 5), txc(1, 6), hq(1, 3), mean)
                end do

                if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
                    ! Buoyancy as next scalar, current value of counter is=inb_scal_array+1
                    if (stats_buoyancy) then
                        if (buoyancy%type == EQNS_EXPLICIT) then
                            call THERMO_ANELASTIC_BUOYANCY(imax, jmax, kmax, s, hq(1, 1))
                        else
                            wrk1d(1:jmax, 1) = 0.0_wp
                            call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, hq(1, 1), wrk1d)
                        end if
                        dummy = 1.0_wp/froude
                        hq(1:isize_field, 1) = hq(1:isize_field, 1)*dummy

                        hq(1:isize_field, 3) = txc(1:isize_field, 3) ! Pass the pressure
                        call AVG_SCAL_XZ(is, q, s, hq(1, 1), &
                                         txc(1, 1), txc(1, 2), txc(1, 4), txc(1, 5), txc(1, 6), hq(1, 3), mean)

                    end if

                    if (imixture == MIXT_TYPE_AIRWATER) then
                        is = is + 1
                        call THERMO_ANELASTIC_THETA_L(imax, jmax, kmax, s, hq(1, 1))

                        hq(1:isize_field, 3) = txc(1:isize_field, 3) ! Pass the pressure
                        call AVG_SCAL_XZ(is, q, s, hq(1, 1), &
                                         txc(1, 1), txc(1, 2), txc(1, 4), txc(1, 5), txc(1, 6), hq(1, 3), mean)
                    end if
                end if

            end if

            call AVG_FLOW_XZ(q, s, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), hq(1, 1), hq(1, 2), hq(1, 3), &
                             mean)

            ! Lagrange Liquid and Liquid without diffusion
            if (part%type == PART_TYPE_BIL_CLOUD_3 .or. part%type == PART_TYPE_BIL_CLOUD_4) then
                l_txc(:, 1) = 1.0_wp; ! We want density
                call PARTICLE_TO_FIELD(l_q, l_txc, txc(1, 5), wrk3d)

                hq(:, 1) = hq(:, 1) + small_wp
                idummy = inb_part - 3 ! # scalar properties solved in the lagrangian
                do is = inb_scal_array + 1 + 1, inb_scal_array + 1 + idummy
                    schmidt(is) = schmidt(1)
                    call PARTICLE_TO_FIELD(l_q, l_q(1, 3 + is - inb_scal_array - 1), hq(1, 2), wrk3d)
                    hq(:, 2) = hq(:, 2)/hq(:, 1)
                    call AVG_SCAL_XZ(is, q, s, hq(1, 2), &
                                     txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), mean)
                end do
            end if

            if (part%type /= PART_TYPE_NONE .and. particle_pdf_calc) then                ! Save particle pathlines for particle_pdf
                write (fname, *) itime; fname = "particle_pdf."//trim(adjustl(fname))
                call PARTICLE_PDF(fname, s, l_g, l_q, l_txc)
            end if

            if (part%type == PART_TYPE_BIL_CLOUD_4) then  ! Save particle residence times
                write (fname, *) itime; fname = "residence_pdf."//trim(adjustl(fname))
                call PARTICLE_RESIDENCE_PDF(fname, l_g%np, l_q)
            end if

        end if

#ifdef TRACE_ON
        call TLAB_WRITE_ASCII(tfile, 'LEAVING STATS_TEMPORAL_LAYER')
#endif

        return
    end subroutine STATISTICS_TEMPORAL

    ! ###################################################################
    ! ###################################################################
    subroutine STATISTICS_SPATIAL()

#ifdef TRACE_ON
        use TLAB_CONSTANTS, only: tfile
        use TLAB_PROCS, only: TLAB_WRITE_ASCII
#endif
        use TLAB_VARS
        use TLAB_ARRAYS
        use DNS_LOCAL
        use BOUNDARY_BUFFER
#ifdef USE_MPI
        use TLAB_MPI_VARS
#endif

        implicit none

        ! -----------------------------------------------------------------------
        integer(wi) is, buff_u_jmin, buff_u_jmax, isize_txc

        ! #######################################################################
#ifdef TRACE_ON
        call TLAB_WRITE_ASCII(tfile, 'ENTERING STATS_SPATIAL_LAYER')
#endif

        ! #######################################################################
        ! Averages
        ! #######################################################################
        if (stats_averages) then
#ifdef USE_MPI
            if (ims_pro == 0) then
#endif
                isize_txc = inb_txc*isize_txc_field

                buff_u_jmin = BuffFlowJmax%size
                buff_u_jmax = jmax - BuffFlowJmax%size + 1
                call AVG_FLOW_SPATIAL_LAYER(isize_txc, buff_u_jmin, buff_u_jmax, mean_flow, txc)

                if (scal_on) then
                    do is = 1, inb_scal
                        call AVG_SCAL_SPATIAL_LAYER(is, isize_txc, buff_u_jmin, buff_u_jmax, mean_flow, mean_scal(1, 1, 1, is), txc)
                    end do
                end if

#ifdef USE_MPI
            end if
#endif
        end if

#ifdef TRACE_ON
        call TLAB_WRITE_ASCII(tfile, 'LEAVING STATS_SPATIAL_LAYER')
#endif

        return
    end subroutine STATISTICS_SPATIAL

end module STATISTICS
