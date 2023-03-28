#include "types.h"
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

    use TLAB_VARS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_PROCS
#endif
    use IO_FIELDS

    implicit none

! -------------------------------------------------------------------
! Grid and associated arrays
    TREAL, dimension(:, :), allocatable, save, target :: x, y, z

! Flow variables
    TREAL, dimension(:, :), allocatable :: q

    TREAL, dimension(:), allocatable :: z1
    TREAL, dimension(:), allocatable :: tmp1, tmp2, tmp3, tmp4, tmp5
    TREAL, dimension(:), allocatable :: wrk1d, wrk2d, wrk3d
    TREAL, dimension(:), allocatable :: profiles

    target q

    TINTEGER ilog
    character*32 fname

    TINTEGER itime_size_max, itime_size, i
    parameter(itime_size_max=128)
    TINTEGER itime_vec(itime_size_max)
    TINTEGER iopt_size_max, iopt_size
    parameter(iopt_size_max=10)
    TREAL opt_vec(iopt_size_max)
    character*512 sRes
    character*32 line
#ifdef USE_MPI
    integer icount
#endif

! Pointers to existing allocated space
    TREAL, dimension(:), pointer :: u, v, w

    TREAL, dimension(:, :), pointer :: dx, dy, dz

! ###################################################################
    call DNS_START

    call IO_READ_GLOBAL('dns.ini')
    call THERMO_INITIALIZE()

#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#endif

    isize_wrk3d = imax*jmax*kmax

! -------------------------------------------------------------------
! allocation of memory space
! -------------------------------------------------------------------
    allocate (x(g(1)%size, g(1)%inb_grid))
    allocate (y(g(2)%size, g(2)%inb_grid))
    allocate (z(g(3)%size, g(3)%inb_grid))

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
        call SCANINICHAR &
            (lfile, 'dns.ini', 'PostProcessing', 'Files', '-1', sRes)
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
        call SCANINICHAR(lfile, 'dns.ini', 'PostProcessing', 'ParamSlCorr', '-1', sRes)
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
    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

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
        call IO_READ_FIELDS(fname, IO_FLOW, imax, jmax, kmax, 3, 0, q)

        write (fname, *) itime; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
        call IO_READ_FIELDS(fname, IO_SCAL, imax, jmax, kmax, inb_scal, inb_scal, z1)

! do correlations
        call SL_CORRELATION_1(ilog, y, dx, dy, dz, u, v, w, z1, profiles, &
                              tmp1, tmp2, tmp3, tmp4, tmp5, wrk1d, wrk2d, wrk3d)

    end do

    call TLAB_STOP(0)
end program SL_CORRELATION
