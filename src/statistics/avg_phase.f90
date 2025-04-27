#include "types.h"

#include "dns_error.h"
#include "dns_const.h"

module AVG_PHASE

    use TLab_WorkFlow
    use TLab_Constants, only: wp, wi, longi, efile
    use TLAB_CONSTANTS, only: sizeofint, sizeofreal
    use FDM, only: g
    use TLab_Memory, only: imax, jmax, kmax, isize_field
    use TLab_Time, only: rtime
    use NavierStokes, only: visc, froude, rossby, prandtl
    use Thermodynamics, only: mach
    use NavierStokes, only: nse_eqns, DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL
    use TLab_Memory, only: inb_flow, inb_scal
    use TLAB_ARRAYS, only: q, s
    use TLab_Arrays, only: wrk2d, wrk3d
    use Thermodynamics, only: gama0
    use TLab_Memory, only: Tlab_Allocate_Real_LONG
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc

    implicit none
    type phaseavg_dt
        sequence
        logical :: active
        integer(wi) :: stride
        character(32) :: type
    end type phaseavg_dt

    interface AvgPhaseSpace
        module procedure AvgPhaseSpaceFieldPtr, AvgPhaseSpaceIndex
    end interface AvgPhaseSpace

    type(phaseavg_dt) :: PhAvg
    real(wp), dimension(:), allocatable, target :: avg_flow, avg_stress, avg_p, avg_scal
    integer(wi) :: nxy, nxz, nyz, nz_total
    integer(wi) :: avg_planes
    character(len=32), parameter :: avgu_name = 'avg_flow'
    character(len=32), parameter :: avgstr_name = 'avg_stress'
    character(len=32), parameter :: avgp_name = 'avg_p'
    character(len=32), parameter :: avgs_name = 'avg_scal'

    integer, parameter, public :: IO_SCAL = 1       ! Header of scalar field
    integer, parameter, public :: IO_FLOW = 2       ! Header of flow field

    public :: AvgPhaseSpace
    public :: avg_flow, avg_p, avg_scal, avg_stress, avg_planes
    public :: PhAvg
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine AvgPhaseInitializeMemory(C_FILE_LOC, restart)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_MPI
        use mpi_f08
        use TLabMPI_VARS, only: ims_pro, ims_npro_i, ims_pro_k
#endif

        implicit none
        character(len=*), intent(in) :: C_FILE_LOC
        integer(wi), intent(in) :: restart
        integer(longi) :: alloc_size
        ! ================================================================== !

        nxy = imax*jmax
        nyz = jmax*kmax
        nxz = imax*kmax

        if (restart == -1) then ! used for calls from outside dns_main / dns.x
            avg_planes = 0
        elseif (mod(restart, PhAvg%stride) == 0) then
            avg_planes = restart/PhAvg%stride
        else
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Number of average planes not an integer. Change stride.')
            call TLAB_STOP(DNS_ERROR_AVG_PHASE)
        end if

#ifdef USE_MPI
        if (ims_pro_k == 0) then
#endif
            alloc_size = imax*jmax*(avg_planes + 1)
            call Tlab_Allocate_Real_LONG(C_FILE_LOC, avg_flow, [alloc_size*inb_flow], 'avgflow.')
            call Tlab_Allocate_Real_LONG(C_FILE_LOC, avg_stress, [alloc_size*6], 'avgstr.') ! allocated not yet coded
            call Tlab_Allocate_Real_LONG(C_FILE_LOC, avg_p, [alloc_size*1], 'avgp.')
            call Tlab_Allocate_Real_LONG(C_FILE_LOC, avg_scal, [alloc_size*inb_scal], 'avgscal.')

            avg_flow(:) = 0.0_wp
            avg_stress(:) = 0.0_wp
            avg_p(:) = 0.0_wp
            avg_scal(:) = 0.0_wp

#ifdef USE_MPI
        end if
#endif
        return
    end subroutine AvgPhaseInitializeMemory

    subroutine AvgPhaseSpaceFieldPtr(localsum, nfield, itr, it_first, it_save, field)
        implicit none
        real(wp), dimension(imax, jmax), intent(inout) :: localsum
        integer(wi), intent(in) :: nfield
        integer(wi), intent(in) :: itr, it_first, it_save
        real(wp), pointer, intent(in) :: field(:)

        integer :: index_loc = 4 ! Index needs to be set to appropriate value for pressure. Needed later for if statement

        call AvgPhaseSpaceExec(localsum, nfield, itr, it_first, it_save, index_loc, field)
    end subroutine AvgPhaseSpaceFieldPtr

    subroutine AvgPhaseSpaceIndex(localsum, nfield, itr, it_first, it_save, index)
        implicit none
        real(wp), dimension(imax, jmax), intent(inout) :: localsum
        integer(wi), intent(in) :: nfield
        integer(wi), intent(in) :: itr, it_first, it_save, index
        real(wp), pointer :: field_loc(:) => null()

        call AvgPhaseSpaceExec(localsum, nfield, itr, it_first, it_save, index, field_loc)
    end subroutine AvgPhaseSpaceIndex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine AvgPhaseSpaceExec(localsum, nfield, itr, it_first, it_save, index, field)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_MPI
        use mpi_f08
        use TLabMPI_VARS, only: ims_comm_z, ims_err, ims_pro, ims_pro_k
#endif

        implicit none
        real(wp), dimension(imax*jmax), intent(inout) :: localsum
        integer(wi), intent(in) :: nfield
        integer(wi), intent(in) :: itr, it_first, it_save, index
        real(wp), dimension(imax, jmax, kmax, 1), target, intent(in) :: field

        integer(wi) :: k, ifld, plane_id
        real(wp), dimension(:), pointer :: avg_ptr
        real(wp), dimension(:), pointer :: loc_field
        integer(wi) :: ipl_srt, ipl_end, iavg_srt, iavg_end, lpl_srt, lpl_end
        ! ================================================================== !
        ! Calculation of the plane id to write the spatial average
        plane_id = 1
        if (it_save /= 0) plane_id = mod((itr - 1) - (it_first), it_save) + 1

        ! Determing the tendency to be written
        if (index == 1) then
            avg_ptr => avg_flow
            call c_f_pointer(c_loc(q), loc_field, shape=[imax*jmax*kmax*nfield])
        elseif (index == 2) then
            avg_ptr => avg_scal
            call c_f_pointer(c_loc(s), loc_field, shape=[imax*jmax*kmax*nfield])
        elseif (index == 4) then
            avg_ptr => avg_p
            call c_f_pointer(c_loc(field), loc_field, shape=[imax*jmax*kmax*nfield])
        elseif (index == 5) then
            avg_ptr => avg_stress
            call c_f_pointer(c_loc(q), loc_field, shape=[imax*jmax*kmax*nfield])
            ! Not yet coded
        else
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Unassigned case type check the index of the field in AvgPhaseSpaceExec')
            call TLAB_STOP(DNS_ERROR_AVG_PHASE)
        end if

        if ((index == 1) .or. (index == 2) .or. (index == 4)) then
            do ifld = 1, nfield
                localsum = 0.0_wp
                ! Computing the space average
                do k = 1, kmax
                    ! Computing the start and end of the plane in field
                    ipl_srt = (ifld - 1)*isize_field + nxy*(k - 1) + 1
                    ipl_end = (ifld - 1)*isize_field + nxy*k

                    localsum = localsum + loc_field(ipl_srt:ipl_end)/g(3)%size !loc_field(:,:,k,ifld)/g(3)%size
                end do

                ! Computing the local sum from start and end of the field for accumulating the space averages
                iavg_srt = (ifld - 1)*nxy*(avg_planes + 1) + nxy*(plane_id - 1) + 1
                iavg_end = (ifld - 1)*nxy*(avg_planes + 1) + nxy*plane_id
#ifdef USE_MPI
                if (ims_pro_k == 0) then
                    call MPI_Reduce(localsum, avg_ptr(iavg_srt:iavg_end), nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err) ! avg_ptr(imax*jmax*restarts*fld)
                else
                    call MPI_Reduce(localsum, MPI_IN_PLACE, nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
                end if
#else
                avg_ptr(iavg_srt:iavg_end) = localsum
#endif

                lpl_srt = (ifld - 1)*nxy*(avg_planes + 1) + nxy*avg_planes + 1
                lpl_end = ifld*nxy*(avg_planes + 1)

#ifdef USE_MPI
                if (ims_pro_k == 0) then
#endif
                    avg_ptr(lpl_srt:lpl_end) = avg_ptr(lpl_srt:lpl_end) + avg_ptr(iavg_srt:iavg_end)/avg_planes
#ifdef USE_MPI
                end if
#endif
            end do
        end if
        return
    end subroutine AvgPhaseSpaceExec

    subroutine AvgPhaseStress(q, itr, it_first, it_save)
        real(wp), dimension(:, :), intent(in) :: q
        integer(wi), intent(in) :: itr
        integer(wi), intent(in) :: it_first
        integer(wi), intent(in) :: it_save

        real(wp), dimension(:), pointer :: u, v, w
        integer(wi) :: plane_id

        target q

        u => q(:, 1)
        v => q(:, 2)
        w => q(:, 3)

        ! Order of computation to reuse cache
        ! uu 1
        ! uv 2
        ! uw 5
        ! vv 3
        ! vw 4
        ! ww 6

        plane_id = 1
        if (it_save /= 0) plane_id = mod((itr - 1) - (it_first), it_save) + 1

        call AvgPhaseCalcStress(u, u, 1, plane_id)
        call AvgPhaseCalcStress(u, v, 2, plane_id)
        call AvgPhaseCalcStress(v, v, 4, plane_id)
        call AvgPhaseCalcStress(v, w, 5, plane_id)
        call AvgPhaseCalcStress(u, w, 3, plane_id)
        call AvgPhaseCalcStress(w, w, 6, plane_id)

    end subroutine AvgPhaseStress

    subroutine AvgPhaseCalcStress(field1, field2, stress_id, plane_id)
#ifdef USE_MPI
        use mpi_f08
        use TLabMPI_VARS, only: ims_comm_z, ims_err, ims_pro, ims_pro_k
#endif
        real(wp), pointer, intent(in) :: field1(:)
        real(wp), pointer, intent(in) :: field2(:)
        integer(wi), intent(in) :: stress_id
        integer(wi), intent(in) :: plane_id

        integer(wi) :: k, ipl_srt, ipl_end, iavg_srt, iavg_end, lpl_srt, lpl_end

        wrk3d(:) = 0.0_wp
        wrk2d(:, :) = 0.0_wp

        wrk3d(:) = field1(:)*field2(:)

        do k = 1, kmax
            ipl_srt = nxy*(k - 1) + 1
            ipl_end = nxy*k

            wrk2d(:, 1) = wrk2d(:, 1) + (wrk3d(ipl_srt:ipl_end))/g(3)%size
        end do

        iavg_srt = (stress_id - 1)*nxy*(avg_planes + 1) + nxy*(plane_id - 1) + 1
        iavg_end = (stress_id - 1)*nxy*(avg_planes + 1) + nxy*(plane_id)

#ifdef USE_MPI
        if (ims_pro_k == 0) then
            call MPI_Reduce(wrk2d, avg_stress(iavg_srt:iavg_end), nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err) ! avg_ptr(imax*jmax*restarts*fld)
        else
            call MPI_Reduce(wrk2d, MPI_IN_PLACE, nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
        end if
#else
        avg_stress(iavg_srt:iavg_end) = wrk2d(:, 1)
#endif

        lpl_srt = (stress_id - 1)*nxy*(avg_planes + 1) + nxy*avg_planes + 1
        lpl_end = (stress_id)*nxy*(avg_planes + 1)

#ifdef USE_MPI
        if (ims_pro_k == 0) then
#endif
            avg_stress(lpl_srt:lpl_end) = avg_stress(lpl_srt:lpl_end) + avg_stress(iavg_srt:iavg_end)/avg_planes
#ifdef USE_MPI
        end if
#endif

    end subroutine AvgPhaseCalcStress

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine IO_WRITE_HEADER(unit, isize, nx, ny, nz, nt, params)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer, intent(in) :: unit, isize
        integer(wi), intent(in) :: nx, ny, nz, nt
        real(wp), intent(in) :: params(isize)

        ! -------------------------------------------------------------------
        integer(wi) offset

        !########################################################################
        offset = 5*SIZEOFINT + isize*SIZEOFREAL

        write (unit) offset, nx, ny, nz, nt

        if (isize > 0) then   ! do not write params to file if there are none
            write (unit) params(1:isize)
        end if

        return
    end subroutine IO_WRITE_HEADER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine IO_Write_AvgPhase(avg_planes, nfield, iheader, it_save, stride, basename, index, avg_ptr, avg_start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use TLab_Memory, only: imax, jmax
        use FDM, only: g
        use TLab_Time, only: rtime, itime
        use NavierStokes, only: visc, froude, rossby, prandtl
        use Thermodynamics, only: mach
        use NavierStokes, only: nse_eqns
        use TLAB_CONSTANTS, only: sizeofint, sizeofreal
        use Thermodynamics, only: gama0

#ifdef USE_MPI
        use mpi_f08
        use TLabMPI_VARS, only: ims_comm_x, ims_err, ims_npro_i, ims_pro_i, ims_pro, ims_comm_z, ims_pro_k
#endif
        implicit none
        integer(wi), intent(in) :: avg_planes
        integer(wi), intent(in) :: nfield
        integer(wi), intent(in) :: it_save
        integer(wi), intent(in) :: stride
        character(len=*), intent(in) :: basename
        integer(wi), intent(in) :: index
        real(wp), dimension(:), intent(in) :: avg_ptr
        integer(wi), intent(in), optional :: avg_start

        character(len=128) :: name
        character(len=32) :: varname(1)
        integer(wi), parameter :: isize_max = 20
        real(wp) :: params(isize_max)
        integer(wi) :: isize, iheader, ifld, ifld_srt, ifld_end
        character(len=10) :: start, end, fld_id
        integer(wi) :: arr_planes, header_offset, ioffset_local
        integer(wi) :: nxy

#ifdef USE_MPI
        integer(kind=MPI_OFFSET_KIND) :: f_offset
        type(MPI_File) :: f_handle
        TYPE(MPI_Datatype) :: ftype, mtype
        type(MPI_Status) :: status
#endif
        nxy = imax*jmax

        if (index > 8 .or. index == 3 .or. index == 5 .or. index == 6 .or. index == 7) then
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Unassigned case type check the index of the field in PhaseAvg_Write')
            call TLAB_STOP(DNS_ERROR_AVG_PHASE)
        end if

        nz_total = it_save/stride + 1

        isize = 0
        isize = isize + 1; params(isize) = rtime
        isize = isize + 1; params(isize) = visc ! inverse of reynolds
        if (iheader == IO_SCAL) then
            isize = isize + 1 + 1                     ! prepare space for schmidt and damkohler

        else if (iheader == IO_FLOW) then
            isize = isize + 1; params(isize) = froude
            isize = isize + 1; params(isize) = rossby
            if (nse_eqns == DNS_EQNS_INTERNAL .or. nse_eqns == DNS_EQNS_TOTAL) then
                isize = isize + 1; params(isize) = gama0
                isize = isize + 1; params(isize) = prandtl
                isize = isize + 1; params(isize) = mach
            end if
        end if
        ! INITIALIZATION OF MPI TYPES SHOULD ONLY BE CARRIED OUT ONCE *AND* DURING INITIALIZATION
        header_offset = 5*SIZEOFINT + isize*SIZEOFREAL

        if (present(avg_start)) then
            write (start, '(I10)') (avg_start)
            write (end, '(I10)') avg_start
        else
            write (start, '(I10)') (itime - it_save + 1)
            write (end, '(I10)') itime
        end if

        do ifld = 1, nfield
            ifld_srt = (ifld - 1)*nxy*(avg_planes + 1) + 1
            ifld_end = ifld*nxy*(avg_planes + 1)
            write (fld_id, '(I10)') ifld
            varname(1) = ''
            if (start == end) then ! write single iteration
                name = trim(adjustl(basename))//trim(adjustl(start))//'.'//trim(adjustl(fld_id))
            else ! write multiple iteration including phase average
                name = trim(adjustl(basename))//trim(adjustl(start)) &
                       //'_'//trim(adjustl(end))//'.'//trim(adjustl(fld_id))
            end if

            arr_planes = (jmax*(avg_planes + 1))
            ! Define the array size for planes and file offset
#ifdef USE_MPI
            f_offset = header_offset + ims_pro_i*imax*8
#endif

#ifdef USE_MPI
            if (ims_pro == 0) then
#endif
#define LOC_STATUS "unknown"
#define LOC_UNIT_ID 75
#include "dns_open_file.h"
                call IO_WRITE_HEADER(LOC_UNIT_ID, isize, g(1)%size, g(2)%size, nz_total, itime, params)
                close (LOC_UNIT_ID)
#ifdef USE_MPI
            end if
#endif

#ifdef USE_MPI
            call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
            if (ims_pro_k == 0) then
                ! Create the MPI derived data types for the file view and contiguous blocks
                call MPI_TYPE_VECTOR(arr_planes, imax, imax*ims_npro_i, MPI_REAL8, ftype, ims_err)
                call MPI_TYPE_COMMIT(ftype, ims_err)
                call MPI_TYPE_CONTIGUOUS(imax, MPI_REAL8, mtype, ims_err)
                call MPI_TYPE_COMMIT(mtype, ims_err)

                ! Open the file for writing
                call MPI_FILE_OPEN(ims_comm_x, name, ior(MPI_MODE_CREATE, MPI_MODE_WRONLY), MPI_INFO_NULL, f_handle, ims_err)

                ! Set the file view
                call MPI_File_set_view(f_handle, f_offset, MPI_REAL8, ftype, 'native', MPI_INFO_NULL, ims_err)

                ! Write the data to the file
                call MPI_FILE_WRITE_ALL(f_handle, avg_ptr(ifld_srt), arr_planes, mtype, status, ims_err)

                ! Close the file
                call MPI_FILE_CLOSE(f_handle, ims_err)

                ! Free the MPI derived data types
                call MPI_TYPE_FREE(ftype, ims_err)
                call MPI_TYPE_FREE(mtype, ims_err)
            end if
#else
#include "dns_open_file.h"
            ioffset_local = header_offset + 1
            write (LOC_UNIT_ID, POS=ioffset_local) avg_ptr(ifld_srt:ifld_end)
            close (LOC_UNIT_ID)
#endif
        end do
        return
    end subroutine IO_Write_AvgPhase

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine AvgPhaseResetVariable()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_MPI
        use mpi_f08
        use TLabMPI_VARS, only: ims_comm_x, ims_err, ims_pro, ims_pro_k
#endif

#ifdef USE_MPI
        if (ims_pro_k == 0) then
#endif
            avg_flow(:) = 0.0_wp
            avg_stress(:) = 0.0_wp
            avg_p(:) = 0.0_wp
            avg_scal(:) = 0.0_wp
#ifdef USE_MPI
        end if
#endif
    end subroutine AvgPhaseResetVariable
end module
