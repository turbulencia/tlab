#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# HISTORY
!#
!# 2007/09/04 - J.P. Mellado
!#              Created
!#
!########################################################################
program SL_CORRELATION
    use TLab_Constants, only: wp, wi

    use TLAB_VARS
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Transpose_Initialize
#endif
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use IO_FIELDS

    implicit none

! -------------------------------------------------------------------
! Grid and associated arrays
    real(wp), dimension(:, :), allocatable, save, target :: x, y, z

! Flow variables
    real(wp), dimension(:, :), allocatable :: q

    real(wp), dimension(:), allocatable :: z1
    real(wp), dimension(:), allocatable :: tmp1, tmp2, tmp3, tmp4, tmp5
    real(wp), dimension(:), allocatable :: wrk1d, wrk2d, wrk3d
    real(wp), dimension(:), allocatable :: profiles

    target q

    integer(wi) ilog
    character*32 fname

    integer(wi) itime_size_max, itime_size, i
    parameter(itime_size_max=128)
    integer(wi) itime_vec(itime_size_max)
    integer(wi) iopt_size_max, iopt_size
    parameter(iopt_size_max=10)
    real(wp) opt_vec(iopt_size_max)
    character*512 sRes
    character*32 line
#ifdef USE_MPI
    integer icount
#endif

! Pointers to existing allocated space
    real(wp), dimension(:), pointer :: u, v, w

    real(wp), dimension(:, :), pointer :: dx, dy, dz

    real(wp) params(1)

    ! ###################################################################
    call DNS_START

    call TLab_Initialize_Parameters('tlab.ini')
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Transpose_Initialize(ifile)
#endif

    call NavierStokes_Initialize_Parameters(ifile)
    call Thermodynamics_Initialize_Parameters(ifile)
    call Gravity_Initialize(ifile)
call Rotation_Initialize(ifile)

! -------------------------------------------------------------------
! allocation of memory space
! -------------------------------------------------------------------
    allocate (q(imax*jmax*kmax, 3))
    allocate (z1(imax*jmax*kmax))
    allocate (tmp1(imax*jmax*kmax))
    allocate (tmp2(imax*jmax*kmax))
    allocate (tmp3(imax*jmax*kmax))
    allocate (tmp4(imax*jmax*kmax))
    allocate (tmp5(imax*jmax*kmax))
    allocate (profiles(jmax*10))
    allocate (wrk1d(isize_wrk1d*5))
    allocate (wrk2d(isize_wrk2d))
    allocate (wrk3d(isize_wrk3d))

! -------------------------------------------------------------------
! File names
! -------------------------------------------------------------------
#ifdef USE_MPI
    if (ims_pro == 0) then
#endif
        call ScanFile_Char &
            (lfile, 'tlab.ini', 'PostProcessing', 'Files', '-1', sRes)
        if (sRes == '-1') then
            write (*, *) 'Integral Iterations ?'
            read (*, '(A512)') sRes
        end if
        itime_size = itime_size_max
        call LIST_INTEGER(sRes, itime_size, itime_vec)
#ifdef USE_MPI
    end if
    call MPI_BCAST(itime_size, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
    icount = itime_size
    call MPI_BCAST(itime_vec, icount, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
#endif

! -------------------------------------------------------------------
! Read local options
! -------------------------------------------------------------------
#ifdef USE_MPI
    if (ims_pro == 0) then
#endif
        call ScanFile_Char(lfile, 'tlab.ini', 'PostProcessing', 'ParamSlCorr', '-1', sRes)
        iopt_size = iopt_size_max
        call LIST_REAL(sRes, iopt_size, opt_vec)

        if (sRes == '-1') then
            write (*, *) 'Use Log Normal Variable (y/n)?'
            read (*, '(A1)') line
            if (line(1:1) == 'y') then; ilog = 1
            else; ilog = 0; end if
        else
            ilog = int(opt_vec(1))
        end if

#ifdef USE_MPI
    end if

    call MPI_BCAST(ilog, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
#endif

! -------------------------------------------------------------------
! Further allocation of memory space
! -------------------------------------------------------------------
! NONE

! -------------------------------------------------------------------
! Read the grid
! -------------------------------------------------------------------
    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, wrk1d(:, 1), wrk1d(:, 2), wrk1d(:, 3))
    call FDM_Initialize(x, g(1), wrk1d(:, 1))
    call FDM_Initialize(y, g(2), wrk1d(:, 2))
    call FDM_Initialize(z, g(3), wrk1d(:, 3))

! ###################################################################
! Define pointers
! ###################################################################
    dx => x(:, 2:) ! to be removed
    dy => y(:, 2:)
    dz => z(:, 2:)

    u => q(:, 1)
    v => q(:, 2)
    w => q(:, 3)

! ###################################################################
! Postprocess given list of files
! ###################################################################
    do i = 1, itime_size

        itime = itime_vec(i)

! read data
        write (fname, *) itime; fname = trim(adjustl(tag_flow))//trim(adjustl(fname))
        call IO_Read_Fields(fname, imax, jmax, kmax, itime, 3, 0, q, params)
        rtime = params(1)

        write (fname, *) itime; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
        call IO_Read_Fields(fname, imax, jmax, kmax, itime, inb_scal, inb_scal, z1, params)

! do correlations
        call SL_CORRELATION_1(ilog, y, dx, dy, dz, u, v, w, z1, profiles, &
                              tmp1, tmp2, tmp3, tmp4, tmp5, wrk1d, wrk2d, wrk3d)

    end do

    call TLab_Stop(0)
end program SL_CORRELATION
