#include "dns_const.h"
#include "dns_error.h"

module PLANES
    use TLab_Constants, only: efile, lfile, wp, wi, fmt_r, small_wp
    use TLab_Memory, only: imax, jmax, kmax, inb_scal_array, inb_flow_array, inb_txc
    use FDM, only: g
    use TLab_Time, only: itime, rtime
    use TLab_Arrays, only: q, s, wrk1d, wrk2d, wrk3d, txc
    use IBM_VARS, only: imode_ibm, ibm_partial
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, scal_on
    use Thermodynamics, only: imixture, MIXT_TYPE_AIRWATER
    use Thermo_Anelastic
    use IO_Fields
    use FI_VORTICITY_EQN
    use FI_GRADIENT_EQN

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

    type(planes_dt), public :: iplanes, jplanes, kplanes

    public :: PLANES_READBLOCK, PLANES_INITIALIZE, PLANES_SAVE

    ! -------------------------------------------------------------------
    integer, parameter :: PLANES_NONE = 0
    integer, parameter :: PLANES_FIX = 1
    integer, parameter :: PLANES_CBL = 2
    integer, parameter :: PLANES_LOG = 3 ! flow fields: log of enstrophy
    ! if scal_on: +log of magnitude of scalar gradient

    character*32 varname(1)
    integer(wi) idummy

    real(wp), pointer :: data_i(:, :, :), data_j(:, :, :), data_k(:, :, :)

    type(io_subarray_dt) :: io_subarray_xy, io_subarray_zy, io_subarray_xz

contains

    ! ###################################################################
    ! ###################################################################
    subroutine PLANES_READBLOCK(bakfile, inifile, block, tag, var)
        character(len=*), intent(in) :: bakfile, inifile, block, tag
        type(planes_dt), intent(out) :: var

        character(len=512) sRes

        ! -------------------------------------------------------------------
        call TLab_Write_ASCII(bakfile, '#'//trim(adjustl(tag))//'=<value>')

        var%n = 0

        call ScanFile_Char(bakfile, inifile, block, trim(adjustl(tag))//'Type', 'fix', sRes)
        if (trim(adjustl(sRes)) == 'none') then; var%type = PLANES_NONE
        elseif (trim(adjustl(sRes)) == 'fix') then; var%type = PLANES_FIX
        elseif (trim(adjustl(sRes)) == 'cbl') then; var%type = PLANES_CBL
        elseif (trim(adjustl(sRes)) == 'log') then; var%type = PLANES_LOG
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong Planes.Type.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        call ScanFile_Char(bakfile, inifile, block, tag, 'void', sRes)
        if (trim(adjustl(sRes)) /= 'void') then
            var%n = MAX_SAVEPLANES; call LIST_INTEGER(sRes, var%n, var%nodes)
        end if

        call ScanFile_Char(bakfile, inifile, block, trim(adjustl(tag))//'Values', 'void', sRes)
        if (trim(adjustl(sRes)) /= 'void') then
            var%n = MAX_SAVEPLANES; call LIST_REAL(sRes, var%n, var%values)
        end if

        if (var%n == 0) var%type = PLANES_NONE ! default

        return
    end subroutine PLANES_READBLOCK

    ! ###################################################################
    ! ###################################################################
    subroutine PLANES_INITIALIZE()
#ifdef USE_MPI
        use mpi_f08, only: MPI_REAL4, MPI_COMM_WORLD
        use TLabMPI_VARS, only: ims_comm_x, ims_comm_z, ims_pro_i, ims_pro_k
#endif

        ! -------------------------------------------------------------------
        integer(wi) inb_scal_dummy

        ! ###################################################################
        if (scal_on) then; inb_scal_dummy = inb_scal_array
        else; inb_scal_dummy = 0; end if

        iplanes%size = (inb_flow_array + inb_scal_dummy + 1)*iplanes%n       ! Flow and scal variables, pressure
        jplanes%size = (inb_flow_array + inb_scal_dummy + 1)*jplanes%n
        kplanes%size = (inb_flow_array + inb_scal_dummy + 1)*kplanes%n

        if (imixture == MIXT_TYPE_AIRWATER) jplanes%size = jplanes%size + 2  ! Add LWP and intgral of TWP

        if (iplanes%type == PLANES_LOG) then; iplanes%size = iplanes%size + iplanes%n
            if (scal_on) iplanes%size = iplanes%size + inb_scal_dummy*iplanes%n
        end if
        if (jplanes%type == PLANES_LOG) then; jplanes%size = jplanes%size + jplanes%n
            if (scal_on) jplanes%size = jplanes%size + inb_scal_dummy*jplanes%n
        end if
        if (kplanes%type == PLANES_LOG) then; kplanes%size = kplanes%size + kplanes%n
            if (scal_on) kplanes%size = kplanes%size + inb_scal_dummy*kplanes%n
        end if

        if (iplanes%size > imax) then
            call TLab_Write_ASCII(efile, 'PLANES_INITIALIZE. Array size imax is insufficient.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
        if (jplanes%size > jmax) then
            call TLab_Write_ASCII(efile, 'PLANES_INITIALIZE. Array size jmax is insufficient.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
        if (kplanes%size > kmax) then
            call TLab_Write_ASCII(efile, 'PLANES_INITIALIZE. Array size kmax is insufficient.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
        if (iplanes%type == PLANES_LOG .or. jplanes%type == PLANES_LOG .or. kplanes%type == PLANES_LOG) then
            if (scal_on) then
                if ((inb_scal_array + 4) > inb_txc) then
                    call TLab_Write_ASCII(efile, 'PLANES_INITIALIZE. Not enough memory for log(grad(scal)) [inb_txc too small].')
                    call TLab_Stop(DNS_ERROR_ALLOC)
                end if
            end if
        end if

        ! Check [ijk]planes%nodes
        if (iplanes%type /= PLANES_NONE) then
            if (any(iplanes%nodes(:iplanes%n) < 1) .or. any(iplanes%nodes(:iplanes%n) > g(1)%size)) then
                call TLab_Write_ASCII(efile, 'PLANES_INITIALIZE. Iplane nodes deceed/exeed grid in x-direction.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
        end if
        if (jplanes%type /= PLANES_NONE) then
            if (any(jplanes%nodes(:jplanes%n) < 1) .or. any(jplanes%nodes(:jplanes%n) > g(2)%size)) then
                call TLab_Write_ASCII(efile, 'PLANES_INITIALIZE. Jplane nodes deceed/exeed grid in y-direction.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
        end if
        if (kplanes%type /= PLANES_NONE) then
            if (any(kplanes%nodes(:kplanes%n) < 1) .or. any(kplanes%nodes(:kplanes%n) > g(3)%size)) then
                call TLab_Write_ASCII(efile, 'PLANES_INITIALIZE. Kplane nodes deceed/exeed grid in z-direction.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
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
            io_subarray_xy%offset = 0
            io_subarray_xy%precision = IO_TYPE_SINGLE
#ifdef USE_MPI
            io_subarray_xy%active = .false.  ! defaults
            if (ims_pro_k == (kplanes%nodes(1)/kmax)) io_subarray_xy%active = .true.
            io_subarray_xy%communicator = ims_comm_x
            io_subarray_xy%subarray = IO_Create_Subarray_XOY(imax, jmax*kplanes%size, MPI_REAL4)
#endif
        end if

        if (iplanes%n > 0) then ! Saving full vertical zOy planes; writing only info of PE containing the first plane
            io_subarray_zy%offset = 0
            io_subarray_zy%precision = IO_TYPE_SINGLE
#ifdef USE_MPI
            io_subarray_zy%active = .false.  ! defaults
            if (ims_pro_i == (iplanes%nodes(1)/imax)) io_subarray_zy%active = .true.
            io_subarray_zy%communicator = ims_comm_z
            io_subarray_zy%subarray = IO_Create_Subarray_ZOY(jmax*iplanes%size, kmax, MPI_REAL4)
#endif
        end if

        if (jplanes%n > 0) then ! Saving full blocks xOz planes for prognostic variables
            io_subarray_xz%offset = 0
            io_subarray_xz%precision = IO_TYPE_SINGLE
#ifdef USE_MPI
            io_subarray_xz%active = .true.
            io_subarray_xz%communicator = MPI_COMM_WORLD
            io_subarray_xz%subarray = IO_Create_Subarray_XOZ(imax, jplanes%size, kmax, MPI_REAL4)
#endif

        end if

        return
    end subroutine PLANES_INITIALIZE

    ! ###################################################################
    ! ###################################################################
    subroutine PLANES_SAVE()
        use TLab_Pointers_3D, only: p_wrk2d
        use TLab_Pointers_3D, only: pointers3d_dt
        use Tlab_Background, only: sbg
        use Averages
        use Integration, only: Int_Simpson

        ! -------------------------------------------------------------------
        integer(wi) offset, j, k, iv, nvars
        character*32 fname, str, fmt
        character*250 line1
        type(pointers3d_dt) :: vars(16)
        real(wp) yrescaled(MAX_SAVEPLANES), henc

        ! ###################################################################
        ! general order of variabeles
        ! [u,v,w,{scal1,...},p,{log(entstrophy),log(grad(scal1)),...)}]

        fmt = '('//fmt_r//')'
        write (line1, fmt) rtime
        write (str, *) itime; str = 'at It'//trim(adjustl(str))//' and time '//trim(adjustl(line1))//'.'

        nvars = 0       ! define pointers
        do iv = 1, inb_flow_array
            nvars = nvars + 1; vars(nvars)%field(1:imax, 1:jmax, 1:kmax) => q(1:imax*jmax*kmax, iv)
        end do
        if (scal_on) then
            do iv = 1, inb_scal_array
                nvars = nvars + 1; vars(nvars)%field(1:imax, 1:jmax, 1:kmax) => s(1:imax*jmax*kmax, iv)
            end do
        end if
        call FI_PRESSURE_BOUSSINESQ(q, s, txc(:, 1), txc(:, 2), txc(:, 3), txc(:, 4), DCMP_TOTAL)
        nvars = nvars + 1; vars(nvars)%field(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 1)

        ! -------------------------------------------------------------------
        if (iplanes%type == PLANES_LOG .or. jplanes%type == PLANES_LOG .or. kplanes%type == PLANES_LOG) then
            call FI_VORTICITY(imax, jmax, kmax, q(:, 1), q(:, 2), q(:, 3), txc(:, 3), txc(:, 4), txc(:, 5))
            txc(1:imax*jmax*kmax, 3) = log10(txc(1:imax*jmax*kmax, 3) + small_wp)
            if (imode_ibm == 1) call IBM_BCS_FIELD(txc(1:imax*jmax*kmax, 3))
            nvars = nvars + 1; vars(nvars)%field(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 3)
            if (scal_on) then
                do iv = 1, inb_scal_array
                    call FI_GRADIENT(imax, jmax, kmax, s(:, iv), txc(:, iv + 3), txc(:, iv + 4))
                    txc(1:imax*jmax*kmax, iv + 3) = log10(txc(1:imax*jmax*kmax, iv + 3) + small_wp)
                    if (imode_ibm == 1) call IBM_BCS_FIELD(txc(1:imax*jmax*kmax, iv + 3))
                    nvars = nvars + 1; vars(nvars)%field(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, iv + 3)
                end do
            end if
        end if

        ! -------------------------------------------------------------------
        if (jplanes%type == PLANES_CBL) then    ! Calculate CBL encroachment height to reevaluate planes
            if (sbg(1)%uslope /= 0.0_wp) then
                do j = 1, jmax
                    wrk1d(j, 1) = AVG1V2D(imax, jmax, kmax, j, 1, s(:, 1))
                end do
                wrk1d(1:jmax, 1) = 2.0_wp*(wrk1d(1:jmax, 1)/sbg(1)%uslope - g(2)%nodes(1:jmax))
                henc = sqrt(Int_Simpson(wrk1d(1:jmax, 1), g(2)%nodes(1:jmax)))
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
            call TLab_Write_ASCII(lfile, trim(adjustl(line1))//' '//trim(adjustl(str)))

            offset = 0
            do iv = 1, nvars
                data_k(:, :, 1 + offset:kplanes%n + offset) = vars(iv)%field(:, :, kplanes%nodes(1:kplanes%n))
                offset = offset + kplanes%n
            end do
            write (fname, *) itime; fname = 'planesK.'//trim(adjustl(fname))
            call IO_Write_Subarray(io_subarray_xy, fname, varname, data_k, kplanes%io)

        end if

        ! -------------------------------------------------------------------
        if (jplanes%n > 0) then
            line1 = 'Writing J-planes'
            do iv = 1, jplanes%n
                write (fname, *) jplanes%nodes(iv); line1 = trim(adjustl(line1))//' '//trim(adjustl(fname))//','
            end do
            call TLab_Write_ASCII(lfile, trim(adjustl(line1))//' '//trim(adjustl(str)))

            offset = 0
            do iv = 1, nvars
                data_j(:, 1 + offset:jplanes%n + offset, :) = vars(iv)%field(:, jplanes%nodes(1:jplanes%n), :)
                offset = offset + jplanes%n
            end do
            if (imixture == MIXT_TYPE_AIRWATER) then    ! Add LWP and integral of total water
                call Thermo_Anelastic_LWP(imax, jmax, kmax, g(2), s(:, inb_scal_array), p_wrk2d, wrk1d, wrk3d)
                data_j(:, 1 + offset, :) = p_wrk2d(:, :, 1)
                offset = offset + 1
                call Thermo_Anelastic_LWP(imax, jmax, kmax, g(2), s(:, inb_scal_array - 1), p_wrk2d, wrk1d, wrk3d)
                data_j(:, 1 + offset, :) = p_wrk2d(:, :, 1)
                offset = offset + 1
            end if
            write (fname, *) itime; fname = 'planesJ.'//trim(adjustl(fname))
            call IO_Write_Subarray(io_subarray_xz, fname, varname, data_j, jplanes%io)

        end if

        ! -------------------------------------------------------------------
        if (iplanes%n > 0) then       ! We transpose to make j-lines together in memory
            line1 = 'Writing I-planes'
            do iv = 1, iplanes%n
                write (fname, *) iplanes%nodes(iv); line1 = trim(adjustl(line1))//' '//trim(adjustl(fname))//','
            end do
            call TLab_Write_ASCII(lfile, trim(adjustl(line1))//' '//trim(adjustl(str)))

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
            call IO_Write_Subarray(io_subarray_zy, fname, varname, data_i, iplanes%io)

        end if

        return
    end subroutine PLANES_SAVE

    !########################################################################
    !########################################################################
    subroutine Thermo_Anelastic_LWP(nx, ny, nz, g, ql, lwp, wrk1d, wrk3d)
        use FDM, only: fdm_dt
        use Integration, only: Int_Simpson

        integer(wi), intent(in) :: nx, ny, nz
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: ql(nx*nz, ny)
        real(wp), intent(out) :: lwp(nx, nz)
        real(wp), intent(INOUT) :: wrk1d(ny)
        real(wp), intent(INOUT) :: wrk3d(nx, ny, nz)

        ! -------------------------------------------------------------------
        integer(wi) i, j, k

        ! ###################################################################
        call Thermo_Anelastic_WEIGHT_OUTPLACE(nx, ny, nz, rbackground, ql, wrk3d)

        do k = 1, nz
            do i = 1, nx
                do j = 1, ny
                    wrk1d(j) = wrk3d(i, j, k)
                end do
                lwp(i, k) = Int_Simpson(wrk1d(1:ny), g%nodes(1:ny))
            end do
        end do

        return
    end subroutine Thermo_Anelastic_LWP

end module PLANES
