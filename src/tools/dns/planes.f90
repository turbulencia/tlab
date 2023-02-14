#include "dns_const.h"
#include "dns_error.h"

module PLANES
    use TLAB_CONSTANTS, only: lfile, efile, wp, wi
    use TLAB_VARS, only: imax, jmax, kmax, inb_scal_array, inb_flow_array
    use TLAB_VARS, only: rbackground, g
    use TLAB_VARS, only: itime, rtime
    use TLAB_ARRAYS, only: q, s, wrk1d, wrk2d, wrk3d, txc
    use TLAB_PROCS
    use THERMO_VARS, only: imixture
    use IO_FIELDS
    implicit none
    save
    private

    integer(wi), parameter, public :: MAX_SAVEPLANES = 20

    type planes_dt
        sequence
        integer(wi) n
        integer(wi) nodes(MAX_SAVEPLANES)
        integer(wi) size
        integer(wi) io(5)
    end type planes_dt

    type(planes_dt), public :: iplanes, jplanes, kplanes
    character*32 varname(1)
    integer(wi) idummy

    real(wp), pointer :: data_i(:, :, :), data_j(:, :, :), data_k(:, :, :)

    public :: PLANES_INITIALIZE, PLANES_SAVE

contains

    ! ###################################################################
    ! ###################################################################
    subroutine PLANES_INITIALIZE()
#ifdef USE_MPI
        use MPI
        use TLAB_MPI_VARS
#endif

        ! -----------------------------------------------------------------------
        integer(wi) id
#ifdef USE_MPI
        integer(wi) :: ndims
        integer(wi), dimension(3) :: sizes, locsize, offset
#endif

        ! ###################################################################
        iplanes%size = (inb_flow_array + inb_scal_array + 1)*iplanes%n           ! Flow and scal variables, pressure
        jplanes%size = (inb_flow_array + inb_scal_array + 1)*jplanes%n
        kplanes%size = (inb_flow_array + inb_scal_array + 1)*kplanes%n
        if (imixture == MIXT_TYPE_AIRWATER) jplanes%size = jplanes%size + 2  ! Add LWP and intgral of TWP

        if (iplanes%size > imax) then
            call TLAB_WRITE_ASCII(efile, 'PLANES_INITIALIZE. Array size imax is is insufficient.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
        if (jplanes%size > jmax) then
            call TLAB_WRITE_ASCII(efile, 'PLANES_INITIALIZE. Array size jmax is is insufficient.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
        if (kplanes%size > kmax) then
            call TLAB_WRITE_ASCII(efile, 'PLANES_INITIALIZE. Array size kmax is is insufficient.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

        ! Pointers
        data_i(1:jmax, 1:iplanes%size, 1:kmax) => txc(1:jmax*iplanes%size*kmax, 2)
        data_j(1:imax, 1:jplanes%size, 1:kmax) => txc(1:imax*jplanes%size*kmax, 2)
        data_k(1:imax, 1:jmax, 1:kplanes%size) => txc(1:imax*jmax*kplanes%size, 2)

        ! Info for IO routines: total size, lower bound, upper bound, stride, # variables
        idummy = jmax*kmax*iplanes%size; iplanes%io = [idummy, 1, idummy, 1, 1]
        idummy = kmax*imax*jplanes%size; jplanes%io = [idummy, 1, idummy, 1, 1]
        idummy = imax*jmax*kplanes%size; kplanes%io = [idummy, 1, idummy, 1, 1]
        varname = ['']

        if (kplanes%n > 0) then ! Saving full vertical xOy planes; writing only info of PE containing the first plane
            id = IO_SUBARRAY_PLANES_XOY
            io_aux(id)%offset = 0
            io_aux(id)%precision = IO_TYPE_SINGLE
#ifdef USE_MPI
            io_aux(id)%active = .false.  ! defaults
            if (ims_pro_k == (kplanes%nodes(1)/kmax)) io_aux(id)%active = .true.
            io_aux(id)%communicator = ims_comm_x
            io_aux(id)%subarray = IO_CREATE_SUBARRAY_XOY(imax, jmax*kplanes%size, MPI_REAL4)
#endif
        end if

        if (iplanes%n > 0) then ! Saving full vertical zOy planes; writing only info of PE containing the first plane
            id = IO_SUBARRAY_PLANES_ZOY
            io_aux(id)%offset = 0
            io_aux(id)%precision = IO_TYPE_SINGLE
#ifdef USE_MPI
            io_aux(id)%active = .false.  ! defaults
            if (ims_pro_i == (iplanes%nodes(1)/imax)) io_aux(id)%active = .true.
            io_aux(id)%communicator = ims_comm_z
            io_aux(id)%subarray = IO_CREATE_SUBARRAY_ZOY(jmax*iplanes%size, kmax, MPI_REAL4)
#endif
        end if

        if (jplanes%n > 0) then ! Saving full blocks xOz planes for prognostic variables
            id = IO_SUBARRAY_PLANES_XOZ
            io_aux(id)%offset = 0
            io_aux(id)%precision = IO_TYPE_SINGLE
#ifdef USE_MPI
            io_aux(id)%active = .true.
            io_aux(id)%communicator = MPI_COMM_WORLD
            io_aux(id)%subarray = IO_CREATE_SUBARRAY_XOZ(imax, jplanes%size, kmax, MPI_REAL4)
#endif

        end if

        return
    end subroutine PLANES_INITIALIZE

    ! ###################################################################
    ! ###################################################################
    subroutine PLANES_SAVE()
        use TLAB_POINTERS_3D

        ! -------------------------------------------------------------------
        integer(wi) offset, j, k
        character*32 fname
        character*250 line1

        ! ###################################################################
        write (fname, *) rtime
        write (line1, *) itime; line1 = 'Writing planes at It'//trim(adjustl(line1))//' and time '//trim(adjustl(fname))//'.'
        call TLAB_WRITE_ASCII(lfile, line1)

        call FI_PRESSURE_BOUSSINESQ(q, s, tmp1, tmp2, tmp3, tmp4)

        if (kplanes%n > 0) then
            offset = 0
            do idummy = 1, inb_flow_array
                data_k(:, :, 1 + offset:kplanes%n + offset) = p_q(:, :, kplanes%nodes(1:kplanes%n), idummy)
                offset = offset + kplanes%n
            end do
            do idummy = 1, inb_scal_array
                data_k(:, :, 1 + offset:kplanes%n + offset) = p_s(:, :, kplanes%nodes(1:kplanes%n), idummy)
                offset = offset + kplanes%n
            end do
            data_k(:, :, 1 + offset:kplanes%n + offset) = tmp1(:, :, kplanes%nodes(1:kplanes%n))
            offset = offset + kplanes%n
            write (fname, *) itime; fname = 'planesK.'//trim(adjustl(fname))
            call IO_WRITE_SUBARRAY(io_aux(IO_SUBARRAY_PLANES_XOY), fname, varname, data_k, kplanes%io)
        end if

        if (jplanes%n > 0) then
            offset = 0
            do idummy = 1, inb_flow_array
                data_j(:, 1 + offset:jplanes%n + offset, :) = p_q(:, jplanes%nodes(1:jplanes%n), :, idummy)
                offset = offset + jplanes%n
            end do
            do idummy = 1, inb_scal_array
                data_j(:, 1 + offset:jplanes%n + offset, :) = p_s(:, jplanes%nodes(1:jplanes%n), :, idummy)
                offset = offset + jplanes%n
            end do
            data_j(:, 1 + offset:jplanes%n + offset, :) = tmp1(:, jplanes%nodes(1:jplanes%n), :)
            offset = offset + jplanes%n
            if (imixture == MIXT_TYPE_AIRWATER) then    ! Add LWP and intgral of TWP
                call THERMO_ANELASTIC_LWP(imax, jmax, kmax, g(2), rbackground, p_s(1, 1, 1, inb_scal_array), p_wrk2d, wrk1d, wrk3d)
                data_j(:, 1 + offset, :) = p_wrk2d(:, :, 1)
                offset = offset + 1
             call THERMO_ANELASTIC_LWP(imax, jmax, kmax, g(2), rbackground, p_s(1, 1, 1, inb_scal_array - 1), p_wrk2d, wrk1d, wrk3d)
                data_j(:, 1 + offset, :) = p_wrk2d(:, :, 1)
                offset = offset + 1
            end if
            write (fname, *) itime; fname = 'planesJ.'//trim(adjustl(fname))
            call IO_WRITE_SUBARRAY(io_aux(IO_SUBARRAY_PLANES_XOZ), fname, varname, data_j, jplanes%io)
        end if

        if (iplanes%n > 0) then       ! We transpose to make j-lines together in memory
            offset = 0
            do idummy = 1, inb_flow_array
                do k = 1, kmax
                    do j = 1, jmax
                        data_i(j, 1 + offset:iplanes%n + offset, k) = p_q(iplanes%nodes(1:iplanes%n), j, k, idummy)
                    end do
                end do
                offset = offset + iplanes%n
            end do
            do idummy = 1, inb_scal_array
                do k = 1, kmax
                    do j = 1, jmax
                        data_i(j, 1 + offset:iplanes%n + offset, k) = p_s(iplanes%nodes(1:iplanes%n), j, k, idummy)
                    end do
                end do
                offset = offset + iplanes%n
            end do
            do k = 1, kmax
                do j = 1, jmax
                    data_i(j, 1 + offset:iplanes%n + offset, k) = tmp1(iplanes%nodes(1:iplanes%n), j, k)
                end do
            end do
            offset = offset + iplanes%n
            write (fname, *) itime; fname = 'planesI.'//trim(adjustl(fname))
            call IO_WRITE_SUBARRAY(io_aux(IO_SUBARRAY_PLANES_ZOY), fname, varname, data_i, iplanes%io)
        end if

        return
    end subroutine PLANES_SAVE

end module PLANES
