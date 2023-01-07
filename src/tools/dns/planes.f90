#include "dns_const.h"
#include "dns_error.h"

module PLANES
    use TLAB_CONSTANTS, only: lfile, efile, wp, wi
    use TLAB_VARS, only: imax, jmax, kmax, isize_field, isize_txc_field, inb_scal_array, inb_flow_array
    use TLAB_VARS, only: rbackground, g
    use TLAB_VARS, only: itime, rtime
    use TLAB_VARS, only: io_aux
    use TLAB_PROCS
    use THERMO_VARS, only: imixture
    use DNS_LOCAL
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

    public :: PLANES_INITIALIZE, PLANES_SAVE

contains

    ! ###################################################################
    ! ###################################################################
    subroutine PLANES_INITIALIZE()
#ifdef USE_MPI
        use MPI
        use TLAB_MPI_VARS
        use IO_FIELDS
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

        ! Info for IO routines: total size, lower bound, upper bound, stride, # variables
        idummy = jmax*kmax*iplanes%size; iplanes%io = [idummy, 1, idummy, 1, 1]
        idummy = kmax*imax*jplanes%size; jplanes%io = [idummy, 1, idummy, 1, 1]
        idummy = imax*jmax*kplanes%size; kplanes%io = [idummy, 1, idummy, 1, 1]
        varname = ['']

        if (kplanes%n > 0) then ! Saving full vertical xOy planes; writing only info of PE containing the first plane
            id = IO_SUBARRAY_PLANES_XOY
            io_aux(id)%offset = 0
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
    subroutine PLANES_SAVE(q, s, p, tmp1, tmp2, tmp3, wrk1d, wrk2d, wrk3d)

        real(wp), intent(IN) :: q(imax, jmax, kmax, inb_flow_array)
        real(wp), intent(IN) :: s(imax, jmax, kmax, inb_scal_array)
        real(wp), intent(INOUT) :: p(imax, jmax, kmax)             ! larger arrays for the Poisson solver,
        real(wp), intent(INOUT) :: tmp1(imax, jmax, kplanes%size)  ! but local shapes are used
        real(wp), intent(INOUT) :: tmp2(imax, jplanes%size, kmax)  ! to arrange plane data.
        real(wp), intent(INOUT) :: wrk3d(jmax, iplanes%size, kmax) ! In this array we transpose data
        real(wp), intent(INOUT) :: tmp3(imax, jmax, kmax, 3)
        real(wp), intent(INOUT) :: wrk1d(*)
        real(wp), intent(INOUT) :: wrk2d(imax, kmax, *)

        ! -------------------------------------------------------------------
        integer(wi) offset, j, k
        character*32 fname
        character*250 line1

        ! ###################################################################
        write (fname, *) rtime
        write (line1, *) itime; line1 = 'Writing planes at It'//trim(adjustl(line1))//' and time '//trim(adjustl(fname))//'.'
        call TLAB_WRITE_ASCII(lfile, line1)

        call FI_PRESSURE_BOUSSINESQ(q, s, p, tmp1, tmp2, tmp3, wrk1d, wrk2d, wrk3d)

        if (kplanes%n > 0) then
            offset = 0
            do idummy = 1, inb_flow_array
                tmp1(:, :, 1 + offset:kplanes%n + offset) = q(:, :, kplanes%nodes(1:kplanes%n), idummy)
                offset = offset + kplanes%n
            end do
            do idummy = 1, inb_scal_array
                tmp1(:, :, 1 + offset:kplanes%n + offset) = s(:, :, kplanes%nodes(1:kplanes%n), idummy)
                offset = offset + kplanes%n
            end do
            tmp1(:, :, 1 + offset:kplanes%n + offset) = p(:, :, kplanes%nodes(1:kplanes%n))
            offset = offset + kplanes%n
            write (fname, *) itime; fname = 'planesK.'//trim(adjustl(fname))
            call IO_WRITE_SUBARRAY4(IO_SUBARRAY_PLANES_XOY, fname, varname, tmp1, kplanes%io, wrk3d)
        end if

        if (jplanes%n > 0) then
            offset = 0
            do idummy = 1, inb_flow_array
                tmp2(:, 1 + offset:jplanes%n + offset, :) = q(:, jplanes%nodes(1:jplanes%n), :, idummy)
                offset = offset + jplanes%n
            end do
            do idummy = 1, inb_scal_array
                tmp2(:, 1 + offset:jplanes%n + offset, :) = s(:, jplanes%nodes(1:jplanes%n), :, idummy)
                offset = offset + jplanes%n
            end do
            tmp2(:, 1 + offset:jplanes%n + offset, :) = p(:, jplanes%nodes(1:jplanes%n), :)
            offset = offset + jplanes%n
            if (imixture == MIXT_TYPE_AIRWATER) then    ! Add LWP and intgral of TWP
                call THERMO_ANELASTIC_LWP(imax, jmax, kmax, g(2), rbackground, s(1, 1, 1, inb_scal_array), wrk2d, wrk1d, wrk3d)
                tmp2(:, 1 + offset, :) = wrk2d(:, :, 1)
                offset = offset + 1
                call THERMO_ANELASTIC_LWP(imax, jmax, kmax, g(2), rbackground, s(1, 1, 1, inb_scal_array - 1), wrk2d, wrk1d, wrk3d)
                tmp2(:, 1 + offset, :) = wrk2d(:, :, 1)
                offset = offset + 1
            end if
            write (fname, *) itime; fname = 'planesJ.'//trim(adjustl(fname))
            call IO_WRITE_SUBARRAY4(IO_SUBARRAY_PLANES_XOZ, fname, varname, tmp2, jplanes%io, wrk3d)
        end if

        if (iplanes%n > 0) then       ! We transpose to make j-lines together in memory
            offset = 0
            do idummy = 1, inb_flow_array
                do k = 1, kmax; do j = 1, jmax
                        wrk3d(j, 1 + offset:iplanes%n + offset, k) = q(iplanes%nodes(1:iplanes%n), j, k, idummy)
                    end do; end do
                offset = offset + iplanes%n
            end do
            do idummy = 1, inb_scal_array
                do k = 1, kmax; do j = 1, jmax
                        wrk3d(j, 1 + offset:iplanes%n + offset, k) = s(iplanes%nodes(1:iplanes%n), j, k, idummy)
                    end do; end do
                offset = offset + iplanes%n
            end do
            do k = 1, kmax; do j = 1, jmax
                    wrk3d(j, 1 + offset:iplanes%n + offset, k) = p(iplanes%nodes(1:iplanes%n), j, k)
                end do; end do
            offset = offset + iplanes%n
            write (fname, *) itime; fname = 'planesI.'//trim(adjustl(fname))
            call IO_WRITE_SUBARRAY4(IO_SUBARRAY_PLANES_ZOY, fname, varname, wrk3d, iplanes%io, tmp1)
        end if

        return
    end subroutine PLANES_SAVE

end module PLANES
