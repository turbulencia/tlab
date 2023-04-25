#include "dns_const.h"
#include "dns_error.h"

module PLANES
    use TLAB_CONSTANTS, only: lfile, efile, wp, wi, fmt_r
    use TLAB_VARS, only: imax, jmax, kmax, inb_scal_array, inb_flow_array
    use TLAB_VARS, only: rbackground, g
    use TLAB_VARS, only: itime, rtime
    use TLAB_ARRAYS, only: q, s, wrk1d, wrk2d, wrk3d, txc
    use TLAB_PROCS
    use THERMO_VARS, only: imixture
    use THERMO_ANELASTIC
    use IO_FIELDS
    implicit none
    save
    private

    integer(wi), parameter :: MAX_SAVEPLANES = 20

    type planes_dt
        sequence
        integer type
        integer(wi) n
        integer(wi) nodes(MAX_SAVEPLANES)
        real(wp) values(MAX_SAVEPLANES)
        integer(wi) size
        integer(wi) io(5)
    end type planes_dt

    type(planes_dt) :: iplanes, jplanes, kplanes
    integer, parameter :: PLANES_NONE = 0
    integer, parameter :: PLANES_FIX = 1
    integer, parameter :: PLANES_CBL = 2

    character*32 varname(1)
    integer(wi) idummy

    real(wp), pointer :: data_i(:, :, :), data_j(:, :, :), data_k(:, :, :)

    public :: iplanes, jplanes, kplanes
    public :: PLANES_READBLOCK, PLANES_INITIALIZE, PLANES_SAVE

contains

    ! ###################################################################
    ! ###################################################################
    subroutine PLANES_READBLOCK(bakfile, inifile, block, tag, var)
        character(len=*), intent(in) :: bakfile, inifile, block, tag
        type(planes_dt), intent(out) :: var

        character(len=512) sRes

        ! -------------------------------------------------------------------
        call TLAB_WRITE_ASCII(bakfile, '#'//trim(adjustl(tag))//'=<value>')

        var%type = PLANES_NONE   ! default
        var%n = 0

        call SCANINICHAR(bakfile, inifile, block, trim(adjustl(tag))//'Type', 'fix', sRes)
        if (trim(adjustl(sRes)) == 'none') then; var%type = PLANES_NONE
        elseif (trim(adjustl(sRes)) == 'fix') then; var%type = PLANES_FIX
        elseif (trim(adjustl(sRes)) == 'cbl') then; var%type = PLANES_CBL
        else
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong Planes.Type.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

        call SCANINICHAR(bakfile, inifile, block, tag, 'void', sRes)
        if (trim(adjustl(sRes)) /= 'void') then
            var%n = MAX_SAVEPLANES; call LIST_INTEGER(sRes, var%n, var%nodes)
        end if

        call SCANINICHAR(bakfile, inifile, block, trim(adjustl(tag))//'Values', 'void', sRes)
        if (trim(adjustl(sRes)) /= 'void') then
            var%n = MAX_SAVEPLANES; call LIST_REAL(sRes, var%n, var%values)
        end if

        return
    end subroutine PLANES_READBLOCK

    ! ###################################################################
    ! ###################################################################
    subroutine PLANES_INITIALIZE()
#ifdef USE_MPI
        use MPI
        use TLAB_MPI_VARS
#endif

        ! -------------------------------------------------------------------
        integer(wi) id

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
        use TLAB_TYPES
        use TLAB_POINTERS_3D, only: p_wrk2d
        use TLAB_VARS, only: sbg
        use AVGS

        ! -------------------------------------------------------------------
        integer(wi) offset, j, k, iv, nvars
        character*32 fname, str, fmt
        character*250 line1
        type(pointers3d_dt) :: vars(16)
        real(wp) yrescaled(MAX_SAVEPLANES), henc
        real(wp) SIMPSON_NU

        ! ###################################################################
        fmt = '('//fmt_r//')'
        write (line1, fmt) rtime
        write (str, *) itime; str = 'at It'//trim(adjustl(str))//' and time '//trim(adjustl(line1))//'.'

        nvars = 0       ! define pointers
        do iv = 1, inb_flow_array
            nvars = nvars + 1; vars(nvars)%field(1:imax, 1:jmax, 1:kmax) => q(1:imax*jmax*kmax, iv)
        end do
        do iv = 1, inb_scal_array
            nvars = nvars + 1; vars(nvars)%field(1:imax, 1:jmax, 1:kmax) => s(1:imax*jmax*kmax, iv)
        end do
        call FI_PRESSURE_BOUSSINESQ(q, s, txc(:, 1), txc(:, 2), txc(:, 3), txc(:, 4))
        nvars = nvars + 1; vars(nvars)%field(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 1)

        ! -------------------------------------------------------------------
        if (jplanes%type == PLANES_CBL) then    ! Calculate CBL encroachment height to reevaluate planes
            if (sbg(1)%uslope /= 0.0_wp) then
                do j = 1, jmax
                    wrk1d(j, 1) = AVG1V2D(imax, jmax, kmax, j, 1, s(:, 1))
                end do
                wrk1d(1:jmax, 1) = 2.0_wp*(wrk1d(1:jmax, 1)/sbg(1)%uslope - g(2)%nodes(1:jmax))
                henc = sqrt(SIMPSON_NU(jmax, wrk1d, g(2)%nodes))
            else
                henc = 1.0_wp
            end if
            yrescaled(1:jplanes%n) = jplanes%values(1:jplanes%n)*henc
            call LOCATE_Y(jplanes%n, yrescaled, jplanes%nodes, g(2)%size, g(2)%nodes)
        end if

        ! -------------------------------------------------------------------
        if (kplanes%n > 0) then
            line1 = 'Writing K-planes'
            do iv = 1, kplanes%n
                write (fname, *) kplanes%nodes(iv); line1 = trim(adjustl(line1))//' '//trim(adjustl(fname))//','
            end do
            call TLAB_WRITE_ASCII(lfile, trim(adjustl(line1))//' '//trim(adjustl(str)))

            offset = 0
            do iv = 1, nvars
                data_k(:, :, 1 + offset:kplanes%n + offset) = vars(iv)%field(:, :, kplanes%nodes(1:kplanes%n))
                offset = offset + kplanes%n
            end do
            write (fname, *) itime; fname = 'planesK.'//trim(adjustl(fname))
            call IO_WRITE_SUBARRAY(io_aux(IO_SUBARRAY_PLANES_XOY), fname, varname, data_k, kplanes%io)

        end if

        ! -------------------------------------------------------------------
        if (jplanes%n > 0) then
            line1 = 'Writing J-planes'
            do iv = 1, jplanes%n
                write (fname, *) jplanes%nodes(iv); line1 = trim(adjustl(line1))//' '//trim(adjustl(fname))//','
            end do
            call TLAB_WRITE_ASCII(lfile, trim(adjustl(line1))//' '//trim(adjustl(str)))

            offset = 0
            do iv = 1, nvars
                data_j(:, 1 + offset:jplanes%n + offset, :) = vars(iv)%field(:, jplanes%nodes(1:jplanes%n), :)
                offset = offset + jplanes%n
            end do
            if (imixture == MIXT_TYPE_AIRWATER) then    ! Add LWP and integral of total water
                call THERMO_ANELASTIC_LWP(imax, jmax, kmax, g(2), rbackground, s(:, inb_scal_array), p_wrk2d, wrk1d, wrk3d)
                data_j(:, 1 + offset, :) = p_wrk2d(:, :, 1)
                offset = offset + 1
                call THERMO_ANELASTIC_LWP(imax, jmax, kmax, g(2), rbackground, s(:, inb_scal_array - 1), p_wrk2d, wrk1d, wrk3d)
                data_j(:, 1 + offset, :) = p_wrk2d(:, :, 1)
                offset = offset + 1
            end if
            write (fname, *) itime; fname = 'planesJ.'//trim(adjustl(fname))
            call IO_WRITE_SUBARRAY(io_aux(IO_SUBARRAY_PLANES_XOZ), fname, varname, data_j, jplanes%io)

        end if

        ! -------------------------------------------------------------------
        if (iplanes%n > 0) then       ! We transpose to make j-lines together in memory
            line1 = 'Writing I-planes'
            do iv = 1, iplanes%n
                write (fname, *) iplanes%nodes(iv); line1 = trim(adjustl(line1))//' '//trim(adjustl(fname))//','
            end do
            call TLAB_WRITE_ASCII(lfile, trim(adjustl(line1))//' '//trim(adjustl(str)))

            offset = 0
            do iv = 1, nvars
                do k = 1, kmax
                    do j = 1, jmax
                        data_i(j, 1 + offset:iplanes%n + offset, k) = vars(iv)%field(iplanes%nodes(1:iplanes%n), j, k)
                    end do
                end do
                offset = offset + iplanes%n
            end do
            write (fname, *) itime; fname = 'planesI.'//trim(adjustl(fname))
            call IO_WRITE_SUBARRAY(io_aux(IO_SUBARRAY_PLANES_ZOY), fname, varname, data_i, iplanes%io)

        end if

        return
    end subroutine PLANES_SAVE

end module PLANES
